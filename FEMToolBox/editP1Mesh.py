"""Interactive editor for creating data accepted by the MATLAB P1Mesh class.

The public entry point is :func:`createP1Mesh`.  It opens a Tk window and
blocks until the user presses Done or closes the window.  The returned arrays
have the same layout as P1Mesh.m:

    nodes     -- float array with shape (2, number_of_nodes)
    triangles -- integer array with shape (3, number_of_triangles), 1-based

Run ``python createP1Mesh.py`` to use the editor as a standalone program.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


MAT_FILE_SIGNATURE = "P1MeshEditor"
MAT_FILE_VERSION = 1


def save_mesh_mat(path: str, nodes: np.ndarray, triangles: np.ndarray) -> None:
    """Write a signed P1Mesh Editor MAT-file."""
    from scipy.io import savemat

    normalized_nodes = _as_nodes(nodes)
    normalized_triangles = _as_triangles(triangles, normalized_nodes.shape[1])
    savemat(path, {
        "p1mesh_editor_format": MAT_FILE_SIGNATURE,
        "p1mesh_editor_version": np.array([[MAT_FILE_VERSION]], dtype=np.int32),
        "nodes": normalized_nodes,
        "triangles": normalized_triangles,
    })


def load_mesh_mat(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Read and validate a signed P1Mesh Editor MAT-file."""
    from scipy.io import loadmat

    data = loadmat(path)
    marker = data.get("p1mesh_editor_format")
    if marker is None:
        raise ValueError("This MAT-file was not created by P1Mesh Editor.")
    marker_text = "".join(str(value) for value in np.asarray(marker).ravel())
    if marker_text != MAT_FILE_SIGNATURE:
        raise ValueError("This MAT-file has an unrecognized P1Mesh Editor signature.")
    version = data.get("p1mesh_editor_version")
    if version is None or np.asarray(version).size != 1:
        raise ValueError("The P1Mesh Editor file version is missing or invalid.")
    if int(np.asarray(version).item()) != MAT_FILE_VERSION:
        raise ValueError("This P1Mesh Editor file version is not supported.")
    if "nodes" not in data or "triangles" not in data:
        raise ValueError("The MAT-file does not contain nodes and triangles.")
    nodes = _as_nodes(data["nodes"])
    triangles = _as_triangles(data["triangles"], nodes.shape[1])
    return nodes, triangles


def _as_nodes(values: Optional[Iterable[Iterable[float]]]) -> np.ndarray:
    if values is None:
        return np.empty((2, 0), dtype=float)
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return np.empty((2, 0), dtype=float)
    if array.ndim != 2 or array.shape[0] != 2:
        raise ValueError("nodes must be a 2-by-N array")
    if not np.isfinite(array).all():
        raise ValueError("node coordinates must be finite")
    return array.copy()


def _as_triangles(
    values: Optional[Iterable[Iterable[int]]], number_of_nodes: int
) -> np.ndarray:
    if values is None:
        return np.empty((3, 0), dtype=np.int64)
    raw = np.asarray(values)
    if raw.size == 0:
        return np.empty((3, 0), dtype=np.int64)
    if raw.ndim != 2 or raw.shape[0] != 3:
        raise ValueError("triangles must be a 3-by-M array")
    array = raw.astype(np.int64)
    if not np.array_equal(raw, array):
        raise ValueError("triangle node indices must be integers")
    if array.min() < 1 or array.max() > number_of_nodes:
        raise ValueError("triangle node index is outside the nodes array")
    return array.copy()


@dataclass
class MeshModel:
    """Small GUI-independent mesh model, using zero-based indices internally."""

    points: list[tuple[float, float]] = field(default_factory=list)
    triangles: list[tuple[int, int, int]] = field(default_factory=list)

    @classmethod
    def from_arrays(cls, nodes: np.ndarray, triangles: np.ndarray) -> "MeshModel":
        return cls(
            [tuple(map(float, nodes[:, i])) for i in range(nodes.shape[1])],
            [tuple(int(v) - 1 for v in triangles[:, i]) for i in range(triangles.shape[1])],
        )

    def add_triangle(self, indices: Iterable[int]) -> None:
        a, b, c = tuple(indices)
        tri = self._oriented_triangle(a, b, c)
        key = frozenset(tri)
        if any(frozenset(existing) == key for existing in self.triangles):
            raise ValueError("this triangle already exists")
        self.triangles.append(tri)

    def _oriented_triangle(self, a: int, b: int, c: int) -> tuple[int, int, int]:
        if len({a, b, c}) != 3:
            raise ValueError("a triangle needs three different nodes")
        pa, pb, pc = self.points[a], self.points[b], self.points[c]
        signed_twice_area = ((pb[0] - pa[0]) * (pc[1] - pa[1])
                             - (pb[1] - pa[1]) * (pc[0] - pa[0]))
        scale = max(1.0, *(abs(v) for point in (pa, pb, pc) for v in point))
        if abs(signed_twice_area) <= 1e-12 * scale * scale:
            raise ValueError("the selected nodes are collinear")
        return (a, b, c) if signed_twice_area > 0 else (a, c, b)

    def _midpoint_node(self, a: int, b: int) -> int:
        midpoint = ((self.points[a][0] + self.points[b][0]) / 2,
                    (self.points[a][1] + self.points[b][1]) / 2)
        tolerance = 1e-12 * max(1.0, abs(midpoint[0]), abs(midpoint[1]))
        for index, point in enumerate(self.points):
            if abs(point[0] - midpoint[0]) <= tolerance and abs(point[1] - midpoint[1]) <= tolerance:
                return index
        self.points.append(midpoint)
        return len(self.points) - 1

    def refine_triangle(self, index: int) -> None:
        """Refine one triangle and remove its edge-sharing neighbors."""
        a, b, c = self.triangles[index]
        ab = self._midpoint_node(a, b)
        bc = self._midpoint_node(b, c)
        ca = self._midpoint_node(c, a)
        children = ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))
        children = [self._oriented_triangle(*tri) for tri in children]
        vertices = {a, b, c}
        refined = []
        for i, triangle in enumerate(self.triangles):
            if i == index:
                refined.extend(children)
            elif len(vertices.intersection(triangle)) < 2:
                refined.append(triangle)
        self.triangles = refined

    def refine_all(self) -> None:
        """Refine every current triangle, reusing midpoints on shared edges."""
        original = list(self.triangles)
        edge_midpoints: dict[tuple[int, int], int] = {}

        def midpoint(a: int, b: int) -> int:
            edge = tuple(sorted((a, b)))
            if edge not in edge_midpoints:
                edge_midpoints[edge] = self._midpoint_node(*edge)
            return edge_midpoints[edge]

        refined = []
        for a, b, c in original:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(self._oriented_triangle(*tri) for tri in
                           ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
        self.triangles = refined

    def delete_node(self, index: int) -> None:
        del self.points[index]
        kept = []
        for triangle in self.triangles:
            if index in triangle:
                continue
            kept.append(tuple(v - (v > index) for v in triangle))
        self.triangles = kept

    def arrays(self) -> tuple[np.ndarray, np.ndarray]:
        nodes = (np.asarray(self.points, dtype=float).T
                 if self.points else np.empty((2, 0), dtype=float))
        triangles = (np.asarray(self.triangles, dtype=np.int64).T + 1
                     if self.triangles else np.empty((3, 0), dtype=np.int64))
        return nodes, triangles


class P1MeshEditor:
    NODE_RADIUS = 6
    DOUBLE_CLICK_INTERVAL_MS = 300

    def __init__(self, model: MeshModel, title: str) -> None:
        import tkinter as tk
        from tkinter import ttk

        self.tk = tk
        self.root = tk.Tk()
        self.root.title(title)
        self.root.geometry("1000x700")
        self.root.minsize(700, 450)
        self.model = model
        self.accepted = False
        self.base_title = title
        self.current_path: Optional[str] = None
        self.mode = "select"
        self.selected_nodes: list[int] = []
        self.selected_triangle: Optional[int] = None
        self.drag_node: Optional[int] = None
        self.drag_recorded = False
        self._last_left_press = None
        self.undo_stack: list[tuple[list[tuple[float, float]], list[tuple[int, int, int]]]] = []
        self.pan_anchor: Optional[tuple[float, float]] = None
        self.offset_x, self.offset_y, self.scale = 500.0, 350.0, 70.0
        self.snap = tk.StringVar(value="0")
        self.status = tk.StringVar()

        bar = ttk.Frame(self.root, padding=6)
        bar.pack(side="top", fill="x")
        # Keep each row short enough for the minimum window width.
        file_bar = ttk.Frame(bar)
        file_bar.pack(fill="x", pady=(0, 4))
        edit_bar = ttk.Frame(bar)
        edit_bar.pack(fill="x", pady=(0, 4))
        mesh_bar = ttk.Frame(bar)
        mesh_bar.pack(fill="x")

        ttk.Button(file_bar, text="Open...", command=self.open_file).pack(side="left", padx=2)
        ttk.Button(file_bar, text="Save", command=self.save_file).pack(side="left", padx=2)
        ttk.Button(file_bar, text="Save As...", command=self.save_file_as).pack(side="left", padx=2)
        ttk.Button(file_bar, text="Cancel", command=self.cancel).pack(side="right", padx=2)
        ttk.Button(file_bar, text="Output...", command=self.confirm_done).pack(side="right", padx=2)
        self.mode_buttons = {}
        for label, mode, key in (("Add Nodes", "node", "N"), ("Add Triangles", "triangle", "T"),
                                 ("Select", "select", "S")):
            button = tk.Button(
                edit_bar, text=f"{label} ({key})", padx=10, pady=3,
                command=lambda m=mode: self.set_mode(m)
            )
            button.pack(side="left", padx=2)
            self.mode_buttons[mode] = button
        ttk.Button(edit_bar, text="Delete", command=self.delete_selected).pack(side="left", padx=(12, 2))
        ttk.Button(mesh_bar, text="Refine Globally", command=self.refine_globally).pack(side="left", padx=2)
        ttk.Button(mesh_bar, text="Fit (F)", command=self.fit_view).pack(side="left", padx=2)
        ttk.Label(mesh_bar, text="Snap spacing:").pack(side="left", padx=(12, 4))
        ttk.Entry(mesh_bar, width=7, textvariable=self.snap).pack(side="left", padx=2)

        self.canvas = tk.Canvas(self.root, background="#fbfbfc", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)
        ttk.Label(self.root, textvariable=self.status, anchor="w", padding=5).pack(fill="x")

        self.canvas.bind("<Configure>", lambda _event: self.redraw())
        self.canvas.bind("<Motion>", self.on_motion)
        self.canvas.bind("<ButtonPress-1>", self.on_left_click)
        self.canvas.bind("<B1-Motion>", self.on_left_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_left_release)
        self.canvas.bind("<ButtonPress-3>", self.on_right_click)
        self.canvas.bind("<ButtonPress-2>", self.on_pan_press)
        self.canvas.bind("<B2-Motion>", self.on_pan_drag)
        self.canvas.bind("<MouseWheel>", self.on_wheel)
        self.canvas.bind("<Button-4>", lambda e: self.zoom(e.x, e.y, 1.15))
        self.canvas.bind("<Button-5>", lambda e: self.zoom(e.x, e.y, 1 / 1.15))
        self.root.bind("<Key>", self.on_key)
        self.root.bind("<Control-z>", self.undo)
        self.root.bind("<Control-Z>", self.undo)
        self.root.protocol("WM_DELETE_WINDOW", self.cancel)
        self.update_mode_buttons()
        self.fit_view()

    def set_mode(self, mode: str) -> None:
        self.mode = mode
        self.selected_nodes.clear()
        self.selected_triangle = None
        self.update_mode_buttons()
        self.redraw()

    def update_mode_buttons(self) -> None:
        for mode, button in self.mode_buttons.items():
            if mode == self.mode:
                button.configure(
                    background="#2563eb", foreground="white",
                    activebackground="#1d4ed8", activeforeground="white",
                    relief="sunken"
                )
            else:
                button.configure(
                    background="#f0f0f0", foreground="#202020",
                    activebackground="#e2e8f0", activeforeground="#202020",
                    relief="raised"
                )

    def to_canvas(self, point: tuple[float, float]) -> tuple[float, float]:
        return self.offset_x + point[0] * self.scale, self.offset_y - point[1] * self.scale

    def to_world(self, x: float, y: float) -> tuple[float, float]:
        point = ((x - self.offset_x) / self.scale, (self.offset_y - y) / self.scale)
        try:
            spacing = float(self.snap.get())
        except ValueError:
            spacing = 0.0
        if spacing > 0:
            point = tuple(round(v / spacing) * spacing for v in point)
        return point

    def nearest_node(self, x: float, y: float) -> Optional[int]:
        best, distance2 = None, (self.NODE_RADIUS + 5) ** 2
        for i, point in enumerate(self.model.points):
            px, py = self.to_canvas(point)
            candidate = (px - x) ** 2 + (py - y) ** 2
            if candidate <= distance2:
                best, distance2 = i, candidate
        return best

    @staticmethod
    def _inside_triangle(p, a, b, c) -> bool:
        def cross(u, v, w):
            return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])
        signs = [cross(a, b, p), cross(b, c, p), cross(c, a, p)]
        return all(v >= 0 for v in signs) or all(v <= 0 for v in signs)

    def triangle_at(self, x: float, y: float) -> Optional[int]:
        p = (x, y)
        for i in reversed(range(len(self.model.triangles))):
            a, b, c = (self.to_canvas(self.model.points[v]) for v in self.model.triangles[i])
            if self._inside_triangle(p, a, b, c):
                return i
        return None

    def on_left_click(self, event) -> None:
        """Detect double clicks with an explicit 300 ms interval."""
        previous = self._last_left_press
        self._last_left_press = (event.time, event.x, event.y, self.mode)
        if previous is not None:
            timestamp, x, y, mode = previous
            if (mode == self.mode
                    and 0 <= event.time - timestamp <= self.DOUBLE_CLICK_INTERVAL_MS
                    and abs(event.x - x) <= 5 and abs(event.y - y) <= 5):
                self._last_left_press = None
                self.on_double_left(event)
                return
        self.on_left_press(event)

    def on_left_press(self, event) -> None:
        node = self.nearest_node(event.x, event.y)
        if self.mode == "node":
            if node is None:
                self.record_undo()
                self.model.points.append(self.to_world(event.x, event.y))
                self.selected_nodes = [len(self.model.points) - 1]
            else:
                self.selected_nodes = [node]
                self.drag_node = node
                self.drag_recorded = False
        elif self.mode == "triangle":
            if node is not None and node not in self.selected_nodes:
                self.selected_nodes.append(node)
                if len(self.selected_nodes) == 3:
                    self.record_undo()
                    try:
                        self.model.add_triangle(self.selected_nodes)
                        self.status.set("Triangle added (orientation normalized counterclockwise).")
                    except ValueError as error:
                        self.undo_stack.pop()
                        self.status.set(str(error))
                    self.selected_nodes.clear()
        else:
            self.selected_nodes = [node] if node is not None else []
            self.selected_triangle = None if node is not None else self.triangle_at(event.x, event.y)
        self.redraw()

    def on_left_drag(self, event) -> None:
        self._last_left_press = None
        if self.mode == "node" and self.drag_node is not None:
            if not self.drag_recorded:
                self.record_undo()
                self.drag_recorded = True
            self.model.points[self.drag_node] = self.to_world(event.x, event.y)
            self.redraw()

    def on_left_release(self, _event) -> None:
        self.drag_node = None
        self.drag_recorded = False

    def on_right_click(self, event) -> str:
        node = self.nearest_node(event.x, event.y)
        menu = self.tk.Menu(self.root, tearoff=False)
        if node is not None:
            self.selected_nodes = [node]
            self.selected_triangle = None
            menu.add_command(label="Set coordinates...",
                             command=lambda: self.edit_node_coordinates(node))
            menu.add_command(label="Delete", command=self.delete_selected)
        else:
            triangle = self.triangle_at(event.x, event.y)
            if triangle is None:
                return "break"
            self.selected_nodes.clear()
            self.selected_triangle = triangle
            menu.add_command(label="Refine into 4 triangles",
                             command=lambda: self.refine_selected_triangle(triangle))
        self.redraw()
        try:
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()
        return "break"

    def on_double_left(self, event) -> str:
        """Edit the coordinates of the node under the pointer."""
        node = self.nearest_node(event.x, event.y)
        if node is None:
            return "break"
        self.drag_node = None
        self.selected_nodes = [node]
        self.selected_triangle = None
        self.redraw()
        self.edit_node_coordinates(node)
        return "break"

    def edit_node_coordinates(self, node: int) -> None:
        from tkinter import messagebox, ttk

        dialog = self.tk.Toplevel(self.root)
        dialog.title(f"Edit node {node + 1}")
        dialog.transient(self.root)
        dialog.resizable(False, False)
        dialog.grab_set()

        x_value = self.tk.StringVar(value=f"{self.model.points[node][0]:.16g}")
        y_value = self.tk.StringVar(value=f"{self.model.points[node][1]:.16g}")
        frame = ttk.Frame(dialog, padding=14)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="x coordinate:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        x_entry = ttk.Entry(frame, width=22, textvariable=x_value)
        x_entry.grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(frame, text="y coordinate:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        y_entry = ttk.Entry(frame, width=22, textvariable=y_value)
        y_entry.grid(row=1, column=1, padx=5, pady=5)
        buttons = ttk.Frame(frame)
        buttons.grid(row=2, column=0, columnspan=2, sticky="e", pady=(10, 0))

        def apply() -> None:
            try:
                x, y = float(x_value.get()), float(y_value.get())
                if not np.isfinite((x, y)).all():
                    raise ValueError
            except ValueError:
                messagebox.showerror(
                    "Invalid coordinates", "Both coordinates must be finite numbers.", parent=dialog
                )
                return
            if self.model.points[node] == (x, y):
                dialog.destroy()
                return
            self.record_undo()
            self.model.points[node] = (x, y)
            dialog.destroy()
            self.redraw()

        ttk.Button(buttons, text="Cancel", command=dialog.destroy).pack(side="right", padx=2)
        ttk.Button(buttons, text="Apply", command=apply).pack(side="right", padx=2)
        dialog.bind("<Return>", lambda _event: apply())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        x_entry.focus_set()
        x_entry.selection_range(0, "end")
        dialog.update_idletasks()
        dialog.geometry(
            f"+{self.root.winfo_rootx() + (self.root.winfo_width() - dialog.winfo_width()) // 2}"
            f"+{self.root.winfo_rooty() + (self.root.winfo_height() - dialog.winfo_height()) // 2}"
        )
        self.root.wait_window(dialog)

    def on_motion(self, event) -> None:
        x, y = self.to_world(event.x, event.y)
        self.status.set(f"Mode: {self.mode}    x={x:.6g}, y={y:.6g}    "
                        f"nodes={len(self.model.points)}, triangles={len(self.model.triangles)}")

    def on_pan_press(self, event) -> None:
        self.pan_anchor = (event.x, event.y)

    def on_pan_drag(self, event) -> None:
        if self.pan_anchor:
            self.offset_x += event.x - self.pan_anchor[0]
            self.offset_y += event.y - self.pan_anchor[1]
            self.pan_anchor = (event.x, event.y)
            self.redraw()

    def on_wheel(self, event) -> None:
        self.zoom(event.x, event.y, 1.15 if event.delta > 0 else 1 / 1.15)

    def zoom(self, x: float, y: float, factor: float) -> None:
        old = self.to_world(x, y)
        self.scale = min(1e7, max(1e-5, self.scale * factor))
        self.offset_x = x - old[0] * self.scale
        self.offset_y = y + old[1] * self.scale
        self.redraw()

    def on_key(self, event) -> None:
        key = event.keysym.lower()
        if key in ("delete", "backspace"):
            self.delete_selected()
        elif key in ("n", "t", "s"):
            self.set_mode({"n": "node", "t": "triangle", "s": "select"}[key])
        elif key == "f":
            self.fit_view()
        elif key == "escape":
            self.selected_nodes.clear(); self.selected_triangle = None; self.redraw()
        elif key in ("return", "kp_enter"):
            self.confirm_done()

    def delete_selected(self) -> None:
        if self.selected_nodes:
            self.record_undo()
            for node in sorted(self.selected_nodes, reverse=True):
                self.model.delete_node(node)
        elif self.selected_triangle is not None:
            self.record_undo()
            del self.model.triangles[self.selected_triangle]
        self.selected_nodes.clear()
        self.selected_triangle = None
        self.redraw()

    def refine_selected_triangle(self, triangle: int) -> None:
        if triangle >= len(self.model.triangles):
            return
        self.record_undo()
        self.model.refine_triangle(triangle)
        self.selected_nodes.clear()
        self.selected_triangle = None
        self.status.set("Triangle refined into four congruent triangles; edge-sharing neighbors removed.")
        self.redraw()

    def refine_globally(self) -> None:
        if not self.model.triangles:
            self.status.set("There are no triangles to refine.")
            return
        old_count = len(self.model.triangles)
        self.record_undo()
        self.model.refine_all()
        self.selected_nodes.clear()
        self.selected_triangle = None
        self.status.set(f"Globally refined {old_count} triangles into {4 * old_count} triangles.")
        self.redraw()

    def record_undo(self) -> None:
        snapshot = (list(self.model.points), list(self.model.triangles))
        self.undo_stack.append(snapshot)
        if len(self.undo_stack) > 100:
            del self.undo_stack[0]

    def undo(self, _event=None) -> str:
        if not self.undo_stack:
            self.status.set("Nothing to undo.")
            return "break"
        points, triangles = self.undo_stack.pop()
        self.model.points = points
        self.model.triangles = triangles
        self.selected_nodes.clear()
        self.selected_triangle = None
        self.drag_node = None
        self.drag_recorded = False
        self.status.set("Undid the last mesh change.")
        self.redraw()
        return "break"

    def fit_view(self) -> None:
        width, height = max(self.canvas.winfo_width(), 700), max(self.canvas.winfo_height(), 400)
        if not self.model.points:
            self.offset_x, self.offset_y, self.scale = width / 2, height / 2, 70.0
        else:
            xs, ys = zip(*self.model.points)
            span_x, span_y = max(xs) - min(xs), max(ys) - min(ys)
            self.scale = min((width - 100) / max(span_x, 1e-9),
                             (height - 100) / max(span_y, 1e-9), 250.0)
            self.offset_x = width / 2 - (min(xs) + max(xs)) * self.scale / 2
            self.offset_y = height / 2 + (min(ys) + max(ys)) * self.scale / 2
        self.redraw()

    def redraw(self) -> None:
        canvas = self.canvas
        canvas.delete("all")
        # Axes provide a stable visual reference without imposing a domain.
        ox, oy = self.to_canvas((0.0, 0.0))
        canvas.create_line(0, oy, canvas.winfo_width(), oy, fill="#e1e4e8")
        canvas.create_line(ox, 0, ox, canvas.winfo_height(), fill="#e1e4e8")
        for i, triangle in enumerate(self.model.triangles):
            coords = [value for node in triangle for value in self.to_canvas(self.model.points[node])]
            fill = "#f9c74f" if i == self.selected_triangle else "#bde0fe"
            canvas.create_polygon(*coords, fill=fill, outline="#3a6ea5", width=2)
        for i, point in enumerate(self.model.points):
            x, y = self.to_canvas(point)
            selected = i in self.selected_nodes
            color = "#e63946" if selected else "#1d3557"
            r = self.NODE_RADIUS + (2 if selected else 0)
            canvas.create_oval(x-r, y-r, x+r, y+r, fill=color, outline="white", width=1)
            canvas.create_text(x+9, y-9, text=str(i+1), anchor="sw", fill="#202020")

    def open_file(self) -> None:
        from tkinter import filedialog, messagebox

        path = filedialog.askopenfilename(
            parent=self.root,
            title="Open P1Mesh Editor file",
            filetypes=(("P1Mesh Editor MAT-files", "*.mat"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            nodes, triangles = load_mesh_mat(path)
        except (OSError, ValueError, TypeError, NotImplementedError) as error:
            messagebox.showerror("Cannot open mesh", str(error), parent=self.root)
            return
        self.model = MeshModel.from_arrays(nodes, triangles)
        self.current_path = path
        self.undo_stack.clear()
        self.selected_nodes.clear()
        self.selected_triangle = None
        self.root.title(f"{self.base_title} — {Path(path).name}")
        self.status.set(f"Opened {path}")
        self.fit_view()

    def _write_file(self, path: str) -> bool:
        from tkinter import messagebox

        try:
            save_mesh_mat(path, *self.model.arrays())
        except (OSError, ValueError, TypeError, ImportError) as error:
            messagebox.showerror("Cannot save mesh", str(error), parent=self.root)
            return False
        self.current_path = path
        self.root.title(f"{self.base_title} — {Path(path).name}")
        self.status.set(f"Saved {path}")
        return True

    def save_file(self) -> bool:
        if self.current_path is None:
            return self.save_file_as()
        return self._write_file(self.current_path)

    def save_file_as(self) -> bool:
        from tkinter import filedialog

        path = filedialog.asksaveasfilename(
            parent=self.root,
            title="Save P1Mesh Editor file",
            defaultextension=".mat",
            filetypes=(("P1Mesh Editor MAT-files", "*.mat"),),
            initialfile=Path(self.current_path).name if self.current_path else "mesh.mat",
        )
        return bool(path) and self._write_file(path)

    def done(self) -> None:
        self.accepted = True
        self.root.destroy()

    def confirm_done(self) -> None:
        from tkinter import ttk

        dialog = self.tk.Toplevel(self.root)
        dialog.title("Output mesh")
        dialog.transient(self.root)
        dialog.resizable(False, False)
        dialog.grab_set()
        frame = ttk.Frame(dialog, padding=18)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="Do you want to save the mesh before output?",
                  font=("TkDefaultFont", 10, "bold")).pack(anchor="w", pady=(0, 14))

        def save_and_output() -> None:
            # Finishing always offers a destination, even if this mesh was saved before.
            if self.save_file_as():
                dialog.destroy()
                self.done()

        def output_only() -> None:
            dialog.destroy()
            self.done()

        ttk.Button(frame, text="Save and Output...", command=save_and_output).pack(fill="x", pady=3)
        ttk.Button(frame, text="Output Without Saving", command=output_only).pack(fill="x", pady=3)
        ttk.Button(frame, text="Continue Editing", command=dialog.destroy).pack(fill="x", pady=(3, 0))
        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.update_idletasks()
        dialog.geometry(
            f"+{self.root.winfo_rootx() + (self.root.winfo_width() - dialog.winfo_width()) // 2}"
            f"+{self.root.winfo_rooty() + (self.root.winfo_height() - dialog.winfo_height()) // 2}"
        )
        dialog.focus_set()
        self.root.wait_window(dialog)

    def cancel(self) -> None:
        self.accepted = False
        self.root.destroy()

    def run(self) -> bool:
        self.root.mainloop()
        return self.accepted


def editP1Mesh(nodes=None, triangles=None, title="P1 mesh editor"):
    """Open the visual editor and return ``(nodes, triangles)``.

    Optional initial data must already use the P1Mesh layout (2-by-N nodes,
    3-by-M one-based triangles).  The editor can open and save its signed MAT
    files.  Choosing either output option returns the edited NumPy arrays;
    closing the window or clicking Cancel returns the initial arrays unchanged.
    """
    initial_nodes = _as_nodes(nodes)
    initial_triangles = _as_triangles(triangles, initial_nodes.shape[1])
    model = MeshModel.from_arrays(initial_nodes, initial_triangles)
    editor = P1MeshEditor(model, str(title))
    if editor.run():
        return editor.model.arrays()
    return initial_nodes, initial_triangles


if __name__ == "__main__":
    result_nodes, result_triangles = editP1Mesh()
    np.set_printoptions(suppress=True)
    print("nodes =")
    print(result_nodes)
    print("triangles =")
    print(result_triangles)

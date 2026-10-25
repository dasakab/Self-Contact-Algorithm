# Visualisation tools (Algorithm2)

- `export_fem_outputs.py`: converts native OUT blocks to per-increment .vtu and .pvd; optional GIF and PNG frames.
- `export_all.sh`: entry point for exporting from the project root.
- `cylinder_live_dashboard_windows.py`: on-demand local HTTP monitor, reads solver log every request.
- `start_live_monitor.sh`: launch the live monitor at http://localhost:8765/monitor.html .
- `adaptive_dashboard.py` and `residual_dashboard.py`: write monitor.html and iteration residual CSV/SVG from live solver logs.
- `watch_second_ring.py`: records candidate inner/outer axial-ring events from accepted solver increments.
- `make_one_line_progress.py`: one line per accepted increment from solver logs.
- `penetration_audit.py`: optional independent triangle intersection test (VTK, NumPy, SciPy).
- `paraview_cross_section_macro.py`: GUI macro to create an X-Z slice inside ParaView.

Typical commands from the project root:

    bash START_LIVE_DASHBOARD.sh
    bash EXPORT_PARAVIEW_GIF.sh
    bash MAKE_INCREMENT_REPORT.sh

Only complete saved increments are included. The `--include-initial` option
adds a reference undeformed configuration translated by the origin in UOPT1/MESH.
It is not a calculated solution increment.

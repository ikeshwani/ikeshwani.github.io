# Isba Keshwani: research site

Static site for GitHub Pages. No build step.

## Publish
1. Create a repo named `ikeshwani.github.io` (public) on GitHub.
2. Put everything in this folder at the repo root (`index.html`, `assets/`, `.nojekyll`).
3. Settings > Pages > Source: "Deploy from a branch", branch `main`, folder `/ (root)`.
4. The site goes live at https://ikeshwani.github.io in a minute or two.

## Preview locally
`python3 -m http.server` in this folder, then open http://localhost:8000

## Files
- `assets/hc-core.js`: the 2D solver (vorticity-streamfunction, runs in the browser)
- `assets/hc-ui.js`: canvas drawing and controls
- `assets/convection.mp4`, `poster.jpg`: cropped Oceananigans.jl movie

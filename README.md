# NF Dataset Dashboards

Self-contained HTML dashboards summarizing NF dataset collections hosted on Synapse. Each dashboard is built from one or more Synapse dataset tables and outputs a single portable HTML file suitable for embedding in a Synapse wiki.

## Repository structure

Each folder corresponds to one Synapse dataset collection and is named `{synapseId}_{short_label}`:

```
syn4939902_jhu_repository/
  build_<name>_dashboard.py    # queries Synapse, processes data, writes HTML
  <name>_dashboard_template.html  # Chart.js / Plotly layout; receives injected data
  <name>_dashboard.html           # generated output — commit after rebuilding
```

The generated HTML is fully self-contained (JS libraries are inlined) so it can be uploaded directly to Synapse as a wiki attachment or a file entity.

## Adding a new dashboard

1. **Create a folder** named `{synapseId}_{short_label}` at the repo root, where `synapseId` is the Synapse ID of the parent folder or project for the dataset.

2. **Copy and adapt the build script** from an existing dashboard. At minimum, update:
   - `DATASETS` — map of Synapse dataset table IDs to display labels
   - `DS_VERSIONS` — pinned version IDs for reproducibility
   - Color/symbol maps to match the new dataset's categories
   - `TEMPLATE` and `OUT_HTML` path constants

3. **Copy and adapt the HTML template.** Update the page `<title>`, header text, and any chart sections that don't apply to the new dataset.

4. **Run the build script** to generate the output HTML:
   ```bash
   cd syn4939902_jhu_repository   # or your new folder
   pip install synapseclient pandas
   python build_<name>_dashboard.py
   ```
   Synapse credentials are read from `~/.synapseConfig` (or environment variables). See the [Synapse Python client docs](https://python-docs.synapse.org/) for login setup.

5. **Commit both the script and the generated HTML** so the dashboard is viewable directly in GitHub and ready to upload to Synapse.

## Embedding on Synapse

Upload `<name>_dashboard.html` to Synapse (as a file entity or wiki attachment), then embed it in a wiki page using an iframe widget:

```
${preview?entityId=synXXXXXXX}
```

or via the Synapse wiki HTML widget pointing at the file's raw URL.

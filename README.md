# Percoll / OptiPrep Density Media Mixer

A GitHub-ready Streamlit app for calculating RBC density-media recipes from Percoll, IPS/IPS74, PBS, OptiPrep, fibrinogen endpoints, and whole-blood additions.

The app implements the calculation branches gathered from the uploaded MATLAB scripts, Percoll protocol collection, OptiPrep data sheet, and thesis-style protocol notes:

- Fast approximate Percoll protocol
- IPS and IPS74 preparation
- SDM from IPS + PBS
- SDM with whole-blood/plasma correction
- OptiPrep working-solution preparation and data-sheet density lookup
- OptiPrep + fibrinogen endpoint preparation
- Two-endpoint density mixing
- Ternary IPS/Percoll + OptiPrep + PBS mixing
- Density correction by adding one stock
- Low/high gradient endpoint helper
- Tube visualization with stacked component fills and volume marks

## Scientific assumptions

The calculator assumes:

1. Additive volumes.
2. Linear volume-weighted density mixing.
3. Mass is `density × volume`.
4. In the whole-blood mode, the target density is the continuous liquid-phase density after plasma/supernatant addition. RBC volume is excluded from that density balance.

For precise RBC work, verify final pH, osmolarity, and density experimentally.

## Local setup

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Each result page also renders a simple tube schematic showing the current mixture composition as stacked liquid layers with ml tick marks. The tube legend wraps long component names and the SVG can be downloaded for documentation or lab notes.

## Repository layout

```text
.
├── app.py                         # Streamlit user interface
├── requirements.txt               # Python dependencies for deployment
├── media_mixer/
│   ├── __init__.py
│   └── calculations.py            # Pure calculation functions
└── tests/
    └── test_calculations.py       # Basic numerical tests
```

## Deployment on Streamlit Community Cloud

1. Create a GitHub repository.
2. Put `app.py`, `requirements.txt`, and the `media_mixer/` package in the repository root.
3. Push to GitHub.
4. In Streamlit Community Cloud, create a new app from that repository.
5. Select `app.py` as the entrypoint file.
6. Deploy.

The repository includes a `requirements.txt` file so the cloud runtime can install the required Python packages.

## Main calculation modes

### SDM from IPS and PBS

```text
p_IPS = (rho_target - rho_PBS) / (rho_IPS - rho_PBS)
p_PBS = 1 - p_IPS
```

### SDM with whole blood

```text
V_RBC = p_RBC * V_final
V_blood = V_RBC / hematocrit
V_plasma = V_blood - V_RBC
V_liquid = V_final - V_RBC
f_plasma = V_plasma / V_liquid

p_IPS = (rho_target - rho_PBS - f_plasma * (rho_plasma - rho_PBS)) / (rho_IPS - rho_PBS)
p_PBS = 1 - p_IPS - f_plasma
```

### Ternary IPS/Percoll + OptiPrep + PBS

```text
O = (rho_target - P*rho_IPS - (1-P)*rho_PBS) / (rho_OptiPrep - rho_PBS)
B = 1 - P - O
```

where `P` is the fixed IPS/Percoll-containing fraction, `O` is the OptiPrep fraction, and `B` is the PBS fraction.


### Tube visualization implementation

The tube graphic is generated as inline SVG and rendered through `streamlit.components.v1.html()`. This avoids Streamlit/Pillow trying to decode SVG bytes as a raster image on some installations.

## Notes for future extensions

Useful next additions:

- Save/load named stock solutions as JSON.
- Add density-meter correction logs.
- Add osmolality prediction for non-standard buffers.
- Add gradient pump CSV export for syringe-pump or stepper-motor protocols.
- Add multi-endpoint batch generation for many target densities.

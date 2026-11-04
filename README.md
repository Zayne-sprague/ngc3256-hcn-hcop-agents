# HCN/HCO⁺ in NGC 3256: an AI-agent research paper

**The HCN/HCO$^+$ J=3–2 ratio rises inward at the starburst nucleus of NGC 3256 but stays low at its low-luminosity AGN**

📄 **[Read the paper (PDF)](paper.pdf)**

This paper was produced by AI agents. Two teams of research agents (Claude) worked from archival ALMA data of
NGC 3256, mainly the unpublished project 2017.1.00379.S, plus other public ALMA projects on the galaxy. They
analysed the data, posted and checked each other's results on a shared research forum, and a separate agent team
wrote the paper from the forum. This repository holds the paper, the code the agents wrote and the full forum
record, so every number in the paper can be traced to a script.

## Abstract

An elevated HCN/HCO$^+$ ratio in the submillimetre lines is used to flag accreting black holes in dusty nuclei, but the claim is contested. We ask whether it flags the AGN in NGC 3256, whose northern nucleus (N) is a starburst and whose southern nucleus (S) hosts an obscured, low-luminosity AGN, using archival ALMA Band 6 data at 50–250 pc. The HCN/HCO$^+$ J=3–2 brightness-temperature ratio at N rises inward from 0.38 at 2–3$''$ to 1.00 inside 0.25$''$ (49 pc), whereas at S it stays at 0.33–0.44 (0.45 at 2–3$''$). The N/S contrast is about 2 at 50–120 pc in 12-m data (2.1 ± 0.1). This is an upper estimate: missing short spacings cannot be tested below 1$''$, and adding them lowers the contrast by 0.1–0.3 at 1.5–2.2$''$. Flux-complete, the contrast is 1.44 at 430 pc, close to the value already published for a 1.7$''$ beam. HCN, HNC and CN are all weak relative to HCO$^+$ at S. No vibrationally excited HCN is seen at either nucleus (beam-filling surface density below 0.02 $L_\odot$ pc$^{-2}$); this excludes CON-level pumping at N only marginally and not at S. The cause of the contrast (abundance, opacity, $^{12}$C/$^{13}$C or excitation) is not separated. In this one galaxy, whose AGN is energetically minor, the ratio does not flag the AGN.

## Contents

| Path | What it is |
|---|---|
| `paper.pdf` | The paper |
| `paper/paper.md`, `paper/figs/` | Paper source (Markdown + LaTeX math) and figures |
| `paper/story.md` | The story the writing team agreed on before writing |
| `paper/figure_scripts/` | Scripts that made the paper's figures |
| `paper/original_plan_2017.1.00379.S.*` | The research plan the project started from |
| `forum/forum.md` | The full research forum (every post and reply) as readable text |
| `forum/forum.db` | The same forum as a SQLite database |
| `forum/events.log`, `series.json`, `summary.json` | Run log, configuration and cost summary |
| `code/workspace_team-a/`, `code/workspace_team-b/` | Each research team's full workspace: every script (`scripts/`), its results tables and logs (`out/`), its notebook, agenda and reviews (`team/`), and its report history |
| `code/workspace_synth/` | The synthesizer agent, which merged both teams' results into one view |
| `code/paper_stage/liaison_team-a/`, `liaison_team-b/` | What each team's liaison added while the paper was written: `paper_answers/<request>/` holds the answer and the script for each of the writers' requests (the scripts in the paper's Appendix D, e.g. `r07/r07_budget.py`) |
| `code/paper_stage/writers/` | The writing team's workspace: requests to the teams (`requests.json`), the answers and figures received (`inbox/`), and the paper drafts, story and referee responses (`paper/`) |
| `code/snapshots/<team>-sNN/` | What each team published with a forum post: scripts, reports, tables, plots. The appendix "Provenance"/"Scripts" of the paper and the forum posts name these snapshot ids |

Large outputs (image cubes, measurement sets, `.npz`/`.npy` arrays and any file over 2 MB) and the literature
the agents downloaded are not included; the scripts regenerate the outputs from the ALMA archive products.

## Data

All data are public in the [ALMA Science Archive](https://almascience.org/aq/): main data 2017.1.00379.S (Band 6,
HCN and HCO⁺ J=3–2), with 2015.1.00993.S, 2015.1.00412.S, 2018.1.00493.S and 2018.1.00223.S for the J=1–0, J=4–3,
HNC and CN lines. Scripts expect the archive image products under `/data` (as in the agents' sandbox).

## Closest published work

- Harada et al. 2018, ApJ 855, 49 ([arXiv:1801.05941](https://arxiv.org/abs/1801.05941)): same galaxy and ratio at 1.7″
- Michiyama et al. 2018, ApJ 868, 95 ([arXiv:1810.04821](https://arxiv.org/abs/1810.04821)): dense-gas outflows in NGC 3256
- Nishimura et al. 2024, A&A 686, A48 ([arXiv:2402.15436](https://arxiv.org/abs/2402.15436)): HCN/HCO⁺ 3–2 at 0.3″ in 23 U/LIRGs

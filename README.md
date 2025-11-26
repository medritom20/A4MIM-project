# A4MIM Project

This is a repository for a team project in the course **Algorithms of Matrix Iterative Methods** (taught by Stefano Pozza, Ph.D., Faculty of Mathematics and Physics, Charles University).

## Team
- **Siqing** – theory lead and presentation
- **Guilherme** – coding lead, implementation of `Algorithm 4 (HS-BCG)`
- **Kamil** – TBD (experiments / validation)
- **Tomas** – implementation of `Algorithm 8 (PDP-BCG)`
- **Jiri** – implementation of `Algorithm 7 (PDR-BCG)`

## Goal
Replicate and experiment with the methods from:

> P. Tichý, G. Meurant, D. Šimonová, **"Block CG algorithms revisited"**, *Numerical Algorithms*, 2025.  
> DOI: [10.1007/s11075-025-02038-4](https://doi.org/10.1007/s11075-025-02038-4) (OA under the CC BY 4.0 license)

---

## Getting started

### Cloning the repository (with submodules)

1. Open a terminal **in the folder where you want to store the project**.
2. Run **one** of the following commands:

HTTPS:
```bash
git clone --recurse-submodules https://github.com/medritom20/A4MIM-project.git
```

SSH:
```bash
git clone --recurse-submodules git@github.com:medritom20/A4MIM-project.git
```

This will create a new subfolder (by default named `A4MIM-project`).  
Open this folder in your editor/IDE (or `cd` into it if you work from the terminal).

### Repository already cloned (without submodules)

If you cloned the repository earlier without `--recurse-submodules`, run the following **inside the repository folder** to fetch and initialise all submodules:

```bash
git submodule update --init --recursive
```

---

## Presentation

The course presentation is stored in the `presentation/` submodule.

The presentation can also be edited collaboratively in the [Overleaf project](https://www.overleaf.com/1674121169kxqqbczbtscd).  
Anyone with this link can edit the project – please do **not** share it outside the team.

# Contributing to CRISPR-HAWK

Thank you for your interest in contributing to **CRISPR-HAWK**!
We welcome contributions from the community to help improve the codebase, fix bugs, expand documentation, or suggest new features.

Please read the following guidelines carefully to help streamline the process.

---

## Table of Contents

* [Code of Conduct](#code-of-conduct)
* [Getting Started](#getting-started)
* [Types of Contributions](#types-of-contributions)
* [Development Setup](#development-setup)
* [Testing](#testing)
* [Style Guide](#style-guide)
* [Pull Requests](#pull-requests)
* [Reporting Issues](#reporting-issues)
* [Contact](#contact)

---

## 📜 Code of Conduct

By participating in this project, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).
We are committed to fostering a welcoming, respectful, and inclusive environment for all contributors.


## 🚀 Getting Started

To start contributing:

1. **Fork** the repository.
2. **Clone** your fork:

   ```bash
   git clone https://github.com/<your-username>/CRISPR-HAWK.git
   cd CRISPR-HAWK
   ```
3. **Create a branch** for your changes:

   ```bash
   git checkout -b my-feature-branch
   ```

---

## 🧩 Types of Contributions

You can contribute in many ways:

* Report/Fix bugs
* Suggest new features
* Suggest/Improve documentation
* Suggest/Improve tests
* Report reproducible issues


## ⚙️ Development Setup

We recommend using a virtual environment or \[conda/mamba] for isolation.

1. Create a virtual environment:

   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

2. Install the package in editable mode with dev dependencies:

   ```bash
   pip install -e .[dev]
   ```

This installs `pytest`, `black`, and any additional tools needed for testing and linting.


## 🧪 Testing

Before submitting your changes, **please make sure all tests pass**:

```bash
pytest
```

To run a specific test:

```bash
pytest tests/test_utils.py::test_function_name
```

If your contribution includes new functionality, please add corresponding unit tests under the `tests/` directory.

For quick checks, you can also run:

```bash
crisprhawk --help
crisprhawk search --help
```

## 🎨 Style Guide

* Use **PEP 8** as a general style guide.
* Format code with `black` (installed via `[dev]` extras):

  ```bash
  black src/ tests/
  ```
* Use descriptive commit messages.
* Write docstrings for public functions and modules.


## 📥 Pull Requests

When your contribution is ready:

1. Push your branch to your fork:

   ```bash
   git push origin my-feature-branch
   ```

2. Open a **Pull Request (PR)** from your fork to the `main` branch of the official repo.

### Checklist before submitting:

* [ ] All tests pass (`pytest`)
* [ ] Code is formatted (`black`)
* [ ] Relevant unit tests added or updated
* [ ] Documentation updated if necessary
* [ ] PR includes a clear description of the problem and solution

We’ll review your PR and provide feedback as soon as possible. Thank you for your contribution!


## 🐛 Reporting Issues

Found a bug? Have a suggestion?

1. Search the [existing issues](https://github.com/pinellolab/CRISPR-HAWK/issues).
2. If it's new, [open a new issue](https://github.com/pinellolab/CRISPR-HAWK/issues/new) and include:

   * A clear description of the issue or request
   * Steps to reproduce (if applicable)
   * System info: OS, Python version, CRISPR-HAWK version
   * Any error messages or logs (use code blocks)


## 📬 Contact

For any questions, collaboration proposals, or off-topic inquiries, feel free to email the authors:

* Manuel Tognon
  <br>manuel.tognon@univr.it

* Rosalba Giugno
  <br>rosalba.giugno@univr.it

* Luca Pinello
  <br>lpinello@mgh.harvard.edu


Thank you for contributing to **CRISPR-HAWK**! 🧬🔍
Let’s build a better genome editing toolbox — together.

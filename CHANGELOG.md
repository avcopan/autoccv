# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- `autoccv.graph`: molecular graph implementation extracted from `automol.graph`
  (`base` generic NetworkX + Pydantic engine, `mol` molecular graph, `ts`
  transition-state graph).
- `autoccv.ccv`: the CCV reaction-mapping algorithm (`CCV`,
  `all_from_reactants_and_products`).
- `autoccv.element`: periodic-table data extracted from `automatics.element`,
  with `scripts/elements-data.py` to regenerate it via `mendeleev`.

## [0.0.0] - YYYY-MM-DD

### Fixed

- Fix 1
- Fix 2...

### Changed

- Change 1
- Change 2...

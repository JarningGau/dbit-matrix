# Documentation System

This page defines how documentation should be maintained.

## Target Audiences

- [`docs/users/`](../users/): workflow users
- [`docs/developers/`](./): code and contract maintainers
- [`docs/maintenance/`](../maintenance/): regression and release operators
- [`docs/archive/`](../archive/): archived notes and legacy material

## Source Of Truth

- workflow overview and navigation: [`README.md`](../../README.md)
- user run instructions: [`docs/users/`](../users/)
- stage input/output contracts: [`docs/developers/contracts.md`](contracts.md)
- config field definitions: [`docs/developers/config-reference.md`](config-reference.md)
- runner behavior and advanced stage entrypoints: [`docs/developers/runner-reference.md`](runner-reference.md)
- maintainer regression procedures: [`docs/maintenance/`](../maintenance/)
- version history: [`docs/maintenance/changelog.md`](../maintenance/changelog.md)

## Update Rules

When behavior changes:

- update the matching user guide if the user-facing run path changes
- update [`docs/developers/contracts.md`](contracts.md) if inputs, outputs, or stage order change
- update [`docs/developers/config-reference.md`](config-reference.md) if config fields or runner rules change
- update [`docs/developers/runner-reference.md`](runner-reference.md) if command-generation behavior changes
- update maintenance docs if validation steps change
- keep `README.md` aligned with current supported workflows and current milestone language

## Writing Rules

- keep user docs short and task-oriented
- avoid exposing implementation detail to new users unless it changes safe usage
- keep stage names, workflow names, and output paths exact
- avoid duplicating normative facts across multiple documents
- prefer linking to the canonical page instead of copying details

## Legacy Docs

The previous Chinese documentation remains under [`doc/`](../../doc/) during migration. New documentation work should target [`docs/`](../).

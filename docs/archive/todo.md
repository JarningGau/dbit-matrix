# Archived Design Notes

This page preserves internal design notes that are not active user-facing or developer-facing contracts.

## Current Theme

- whether calling should remain CpG-only or expand to CH contexts

## Existing Caller Shape

The current `scripts/methy_caller.py` is a CpG-specific caller rather than a general cytosine caller.

## Design Direction Under Discussion

- do not extend the current CpG dinucleotide logic directly to CH calling
- if CH support is added, prefer a strand-aware cytosine caller
- define the MVP output contract before changing summary or downstream expectations

## Follow-Up Questions

- should the workflow keep a CpG-only contract
- if CH is added, should outputs gain explicit `strand` or `context` columns
- should the caller be refactored from CpG-specific logic to cytosine-event logic

# Archived Progress Notes

This page preserves internal project status notes that do not belong in the main user or maintainer docs.

## Snapshot

- the EMSeq standalone mainline is closed end to end
- the RNA standalone workflow currently covers `demux_extract_bc -> align`
- RNA `align` writes STARsolo matrix outputs under `work/<sample>/solo/`

## Historical Milestones

- M1: `fastp_split + demux_extract_bc`
- M2: `align`
- M3: `pool`
- M4: `split`
- M5-M7: `call`, `mbias`, and `summary`
- later addition: `saturation`

## Open Internal Topics

- evaluate `call` performance and statistics consistency
- decide whether `mbias` should eventually feed trimming or calling defaults
- decide whether `host_mito` should stay a single aggregate output or gain layered reporting

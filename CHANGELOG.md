# phac-nml/genomic_address_service: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).



## v0.2.2 - [2024-08-01]

### `Fixed`

- Removed suffix from created databases, this resolves an issue where databases are not found as the wrong suffix is passed to the database option [PR 36](https://github.com/phac-nml/locidex/pull/36)

## v0.2.1 - [2024-07-29]

### `Fixed`

- Changed default data format for database configs to year-month-day format [PR 32](https://github.com/phac-nml/locidex/pull/32/commits/3afbd16880c9e70738fa73f7adfcedda59825983)

## v0.2.0

### `Added`

- Incorporated both unit tests and end-to-end tests for workflows

- Introduced manifest module for management of multiple databases

### `Fixed`

- Reduced code coupling issue by removing incorporating dataclasses over dictionaries. Changes were made in a backwards compatible fashion.

## v1.0dev - [date]

Initial release of phac-nml/locidex

### `Added`

Added locidex workflow mermaid code to ./assets/

Created CHANGELOG.md

### `Fixed`

Changed README format to standard DAAD README

Changed locidex workflow to  mermaid format

### `Dependencies`

### `Deprecated`

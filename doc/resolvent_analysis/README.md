# Resolvent Analysis Documentation

This directory contains documentation for the resolvent analysis implementation in ADflow.

## Quick Start

See [CRITICAL_SETUP_REQUIREMENT.md](../../CRITICAL_SETUP_REQUIREMENT.md) for setup instructions.

See [IMPLEMENTATION_COMPLETE.md](../../IMPLEMENTATION_COMPLETE.md) for complete overview.

## Implementation Details

- **Theory**: See main implementation files and academic references
- **Usage Examples**: In IMPLEMENTATION_COMPLETE.md
- **Performance**: See SESSION_SUMMARY.md for benchmarks
- **Validation**: All tests in tests/ directory

## Key Files

### Essential (Keep in Root)
- `CRITICAL_SETUP_REQUIREMENT.md` - Must-read setup guide
- `IMPLEMENTATION_COMPLETE.md` - Complete implementation overview
- `SESSION_SUMMARY.md` - Development session summary

### Archived (This Directory)
- Historical development notes
- Intermediate validation reports
- Migration guides (if needed)

## Test Files

All tests are in `tests/` directory:
- `test_resolvent_simple.py` - Explicit method validation
- `test_resolvent_matrix_free.py` - Matrix-free method validation
- `test_adjoint_preconditioner.py` - Preconditioner unit tests
- `test_ilu_preconditioner.py` - ILU preconditioner tests
- `test_resolvent_real_form.py` - Real-valued form tests

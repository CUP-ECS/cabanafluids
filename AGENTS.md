# Working in this repository

Instructions for AI agents working in this repository, and the rules a human owns
when an agent works on their behalf. Reachable as `CLAUDE.md`, a symlink, so
harnesses that look for that name find this file.

This copy is adapted from `COMPASS-PSAAP/compass-template`. The standing rules below
apply to every repository in the organization; change them in the template, not here.
**What this repository holds**, **Repository map**, and **Additional engineering
conventions** describe CabanaFluids.

## Contents

- [What this repository holds](#what-this-repository-holds)
- [Repository map](#repository-map)
- [Repository configuration](#repository-configuration)
- [What agents are for here](#what-agents-are-for-here)
- [Agent instructions](#agent-instructions)
  - [Ask a human first](#ask-a-human-first)
  - [Never](#never)
  - [How to work here](#how-to-work-here)
  - [Git and pull requests](#git-and-pull-requests)
  - [Evidence and reporting](#evidence-and-reporting)
  - [Writing](#writing)
  - [Durable artifacts](#durable-artifacts)
  - [Machines, allocations, and data](#machines-allocations-and-data)
  - [Revision history discipline](#revision-history-discipline)
  - [Memory and session continuity](#memory-and-session-continuity)
- [Human instructions](#human-instructions)
  - [Running agents safely](#running-agents-safely)
  - [The coding feature loop](#the-coding-feature-loop)
  - [What the human owns](#what-the-human-owns)
- [Additional engineering conventions](#additional-engineering-conventions)

## What this repository holds

CabanaFluids is a 2D finite difference solver for incompressible Euler flow that
advects a passive quantity, built on Kokkos and `Cabana::Grid`. It serves as a
benchmark for halo-exchange, neighbor-collective, and GPU communication strategies,
so predictable communication patterns matter as much as numerical accuracy. The
`advection` example in `examples/` drives the library in `src/`. `FUTURE.md` holds
the known bugs, the remaining test plan, and the planned replacement of
semi-Lagrangian advection; consult it before changing the numerics. The repository
lives at `CUP-ECS/CabanaFluids`, and no other repository owns its specification.

## Repository map

| Path | Owns |
|---|---|
| `AGENTS.md` | This file. The agent and the human instructions, and the repository map. Inherited from the template; the repository description, map, and conventions are adaptable, but standing-rule changes go back to the template |
| `CLAUDE.md` | A symlink to `AGENTS.md`. Edit `AGENTS.md`; never replace the symlink with a copy, or the two names drift |
| `README.md` | Human-facing: what this repository is, how to build it, how to use it. Placeholder inherited from the template but changeable |
| `CHANGELOG.md` | Concise history of changes to this repository, newest first. The template's entries are replaced on adoption |
| `LICENSE` | BSD 3-Clause. Changing it is the author's decision, not an agent's. Inherited from the template and frozen |
| `.github/` | Pull-request and issue templates, and the `clang-format` workflow, inherited from the template but changeable. The build-and-test and CI-container workflows are this repository's own. Changing a workflow or adding one is CI configuration and needs a human |
| `.clang-format` | The C and C++ style CI enforces. This repository keeps its own LLVM-based Allman style rather than the template's |
| `.gitignore` | Build, editor, and Spack artifacts kept out of the repository. Inherited from the template, with `build-*` added for Spack build symlinks |
| `spack-envs/llnl/` | Placeholder Spack environments for Dane, Tioga, and Tuolumne. Inherited from the template but changeable |
| `MEMORY.md` | Agent session continuity: current state, decisions the code records only as outcomes, easy-to-miss procedures, and open questions |
| `FUTURE.md` | Known bugs, fixed bugs with their tests, the unit and convergence test plan, and the planned MUSCL and SSP-RK3 advection scheme |
| `CMakeLists.txt`, `cmake/` | The CMake build, the find modules for Silo and `clang-format`, and Cabana's test harness in `cmake/test_harness/` |
| `src/` | The header-mostly `cabanafluids` library: mesh, problem manager, boundary conditions, time integration, pressure projection, and Silo output |
| `tests/` | GoogleTest unit tests, one `tst<Name>.hpp` per suite, each built per Kokkos backend and registered with CTest at 1, 2, and 4 MPI ranks |
| `examples/` | The `advection` driver and its command-line options |
| `configs/` | Spack environments and the CMake host config for CI, generic CUDA systems, class tutorials, and macOS |

## Repository configuration

**New repositories start from `COMPASS-PSAAP/compass-template`** — on
GitHub through "Use this template", or with
`gh repo create COMPASS-PSAAP/<name> --template COMPASS-PSAAP/compass-template`.

**Keep inherited files close to the template.** Drift is then visible as a diff
rather than buried in a rewrite. If possible, repository-specific additions to an 
inherited file go in a clearly marked "For This Repository" block at the end of it.
An improvement that belongs to every repository goes to the template as a pull request 
there.

Parts of the configuration will have nothing to act on. A pure-Python repository has
no use for `.clang-format` or its workflow, and a repository of prose has no use for
the Spack environments. Leaving them in place is the cheaper default; pruning them is
the author's call, and either way the decision is recorded in `CHANGELOG.md`
rather than left to be rediscovered.

## What agents are for here

An agent's output is a **proposal**, and it is not part of this repository until a
human has read it and merged it.

Work agents do well here, unsupervised: reading the repository and answering
questions about it, drafting an implementation on a branch, writing tests against a
stated specification, reproducing a failure, mechanical refactors with a green test
suite on both sides, and drafting prose for a human to argue with.

Work agents do badly here, and must not do alone: deciding what a result means,
choosing a threshold, resolving a disagreement between two documents, deciding that a
failing test was wrong, and anything whose cost lands on someone outside the session
— machine time, a collaborator's inbox, a published page, a merged branch.

## Agent instructions

Directions for AI agents, which the humans working with them should know too.

### Ask a human first

Stop and ask before any of these, and do not treat approval for one as approval for
the next:

- **Pushing to a shared branch**, force-pushing, or rewriting published history
- **Merging a pull request**, approving a review, or closing an issue
- **Submitting a batch job**, reserving nodes, or anything else that spends
  allocation hours
- **Installing into a shared environment** or modifying a shared module or Spack
  installation
- **Adding a dependency**, vendoring third-party code, or changing the license
- **Changing repository settings**, CI configuration, secrets, or branch protection
- **Deleting files, data, or results**, and anything whose undo is a restore from
  backup
- **Sending anything outside the repository** — publishing a page, opening an issue on
  another organization's repository, posting to an external service, or mailing a
  collaborator
- **A load-bearing rewrite.** A correction that invalidates a deliverable, a
  deadline, or a claim changes what the project can say. Write the proposal; apply it
  once it is accepted

### Never

- **Commit a secret** — a token, a key, a password, a certificate, a machine
  hostname-plus-credential pair — or a file that contains one. If a secret has 
  already been committed, stop and tell a human.
- **Send sensitive, export-controlled, proprietary, or unpublished collaborator
  material to an external service without permission.** Unpublished results and
  code shared under an agreement are included. By starting a session the human
  authorizes sending this repository's contents to the model behind it, but that
  does not authorize sending restricted third-party material anywhere. When in
  doubt, treat information as restricted and ask.
- **Invent a result, a number, a date, or a citation.** An unverified reference is
  marked unverified. A number you did not compute is not reported.
- **Make a test pass by weakening it, skipping it, or asserting what the code does.**
  A failing test is a finding.
- **Report work as done that was not run.** If the tests were not run, say so; if a
  step was skipped, say which.
- **Write in someone's voice about what they did.** Text that asserts what the authors
  chose, observed, understood, or learned comes from those people. See
  [Writing](#writing).
- **Edit `LICENSE`**, or replace the `CLAUDE.md` symlink with a copy.

### How to work here

- **Re-read the files at the start of each task.** They change between sessions and
  outside them. Do not act on a remembered version, including one recalled from
  memory; memory records what was true when it was written.
- **Keep user edits.** People edit these files directly. Never revert an edit you did
  not make, and never regenerate a file wholesale when an edit will do.
- **Keep the diff small and reviewable.** A large reflow hides the substantive change
  inside it, and a reviewer who cannot see the change cannot catch the error.
- **Match the surrounding code and prose.** Naming, comment density, wrapping, and
  idiom. A patch that reads as if it came from somewhere else costs a reviewer twice.
- **Never remove a comment**, and never delete a claim to resolve a disagreement.
  Revise it and record the revision. 
- **Retire, do not erase.** When a decision is reversed, say so where the decision
  lives and in `CHANGELOG.md`. Retired decisions that vanish come back.
- **Update the cross-references you invalidate.** Section numbers, filenames, symbol
  names, and issue IDs are cited across files; grep for the identifier you moved.
- **Say what you are unsure of.** An open question written down is cheap; a guess
  hardened into an implementation is not.

### Git and pull requests

- **Get default branch protection enabled.** After the repository is created or
  forked, ask a human to turn it on if the organization has not already turned
  it on by default.
- **Never commit to the default branch.** Branch off `origin/HEAD` with a name that
  says what the branch does, and open a pull request.  This includes initial setup — 
  you may create the default branch (`main`) if it does not exist, but every change 
  to it goes through a pull request.
- **Fill in the pull-request template**, including the checklists. They are the
  organization's review contract.
- **Say in the pull request that an agent wrote it**, and what a reviewer should look
  at hardest.
- **Do not add agent `Co-Authored-By` trailers.** Attribution belongs in the pull
  request description, where it can be read alongside what was reviewed.
- **Prefer amending to stacking fix commits**, and update the message to match what
  the commit now does. Subject under 51 characters, body lines under 73, why over
  what, real newlines.
- **One change per pull request.** A formatting sweep and a behavior change in the
  same branch cannot be reviewed, only skimmed.

### Evidence and reporting

These apply to results, to documents, and to anything an agent asserts in a pull
request.

- **Distinguish measured from inferred, and predicted from observed**, on every
  number. A number derived from a model is labeled with the model and its inputs.
- **Every performance claim carries a trial count and a spread.** Single-run numbers
  are not results.
- **Report the tail, not the mean alone**, wherever the effect could live in the tail,
  and report a dispersion measure alongside every location measure.
- **Say what was dropped.** A truncated list, a sampled set, a capped sweep — a
  bounded analysis presented without its bound reads as full coverage.
- **Record provenance.** What produced the number: the commit, the inputs, the
  machine, the environment, the seed. A number with no provenance block is not a
  result.
- **Do not invent constants.** Machine-specific values default to a published number,
  cited where it is used, and are arguments rather than literals.

### Writing

Prose here is read by people who will argue with it, so it is written to be argued
with. Strunk and White's *The Elements of Style* and Embree's *Gospel*
(https://www.cs.unm.edu/~bridges/embree-gospel.pdf) are the house sources.

- **Orient the reader first, then organize their thoughts.** Assume they will not read
  the whole document or read it in order.
- **Omit needless words, use the active voice, and prefer the specific to the general.**
  Write with nouns and verbs. 
- **Keep the tone of an engineer dressed for work** — matter-of-fact and
  straight-from-the-shoulder, neither stuffy nor cute.
- **Avoid the habits that make prose read as machine-written**: stock openers and
  wrap-up closers, hedge-stacking, empty transitions, reflexive rule-of-three lists,
  bulleting a continuous thought, and "not X, but Y" constructions that bury the point.
- **Do not supply the human's own content.** A paper's contributions, a report on what
  was tried and what was learned, a status update, a response to a reviewer, and a
  progress report all assert what particular people did and concluded. The person
  drafts those; the agent critiques the draft, tightens it, and says what is missing.
  Generated experience is fabricated evidence wearing a different hat.

### Durable artifacts

The session scratchpad is not storage. Anything someone outside the session will read
again — a figure, a generated table, a campaign script, a shareable page — belongs
in this repository, on a branch, pushed. Scratchpad files disappear with the session and
are unrecoverable.

- **A generated artifact is not finished until it is committed.** Publishing puts a
  rendering on a hosted page; committing puts it where someone else can regenerate,
  diff, and correct it. Commit it locally, and ask before publishing.
- **A derived artifact names its sources inside the artifact** — which files, at which
  commit — and states that the sources win where the two disagree. A rendering with no
  provenance becomes a second source of truth that nobody can date.
- **Derived files live in their own directory**, never beside the authoritative ones.
  The directory is the signal that the file is generated rather than owned.
- **One-shot helper scripts are not artifacts.** A patch script already applied, or a
  validator written to check one file once, stays in the scratchpad.

### Machines, allocations, and data

- **Agents do not spend machine time.** Building, running a unit test suite, and
  running a single-node smoke test are ordinary work; submitting to a queue, holding
  nodes, and running a campaign are not, and need a human each time.
- **Where the Spack environments apply, use them** rather than a hand-assembled module
  list, and record any deviation in the pull request.
- **Record the environment with every timing** — machine, allocation, node count,
  compiler, MPI, relevant runtime variables. A timing whose environment is unstated
  cannot be compared to anything.
- **Treat data as read-only** unless the task is to produce it. Do not overwrite a
  results file, and do not delete an input to reclaim space.
- **Cache downloads under a stated directory with a checksum**, and never re-download
  silently.

### Revision history discipline

Each substantive change is recorded concisely in `CHANGELOG.md`, newest first.
Say what changed and why, with the date. Editorial cleanups need no entry.

### Memory and session continuity

A session's memory dies with it, and the next agent inherits the files rather than the
reasoning. When agent reasoning is worth carrying forward but belongs in no file,
write it in a `MEMORY.md` file at the root of the repository. The last step of a
working session is to bring it up to date:

- **Current state** — what was being worked on, and whether it is done, blocked, or
  half finished
- **Design decisions and their reasoning**, where the code or the document records only
  the outcome, and reversals, so they are not proposed again
- **Operational procedures that are easy to get wrong**
- **Open questions and next steps**, including threads waiting on a person or an event
  rather than on work

Keep it short and dated without a revision history, point at the file that owns each 
area rather than duplicating it, and verify before acting on it: it records what was 
true when it was written.  Add its row to the **Repository map** when you create it.

## Human instructions

Directions for the humans running the agents, which the agents also read.

### Running agents safely

- **Keep credentials out of the repository and out of harness configuration.** Inject
  API keys from a password manager or a Codespaces secret at the start of the session.
- **Scope tokens to this repository**, not to the account. A personal access token
  handed to an MCP server has whatever reach you gave it, for as long as the session runs.
- **Add an MCP server only when you can say what it is for.** A GitHub server that
  reads issues, opens pull requests, and checks Actions status earns its place by
  removing hand-copying. Each additional one is another path out of the sandbox.
- **Match the model to the task.** Mechanical work — boilerplate, formatting, a
  straightforward fix, a test case — does not need the expensive model. Watch the context
  too: a session that repeatedly dumps whole files and long command output spends
  budget without adding understanding, so point the agent at line ranges.
- **Assume the next session runs a different model.** Instructions live in `AGENTS.md`,
  which every harness reads, with `CLAUDE.md` as a symlink to it. What a session
  learned goes in the repository, not in one harness's configuration.

### The coding feature loop

One change at a time, and the same way every time:

1. **Branch** off `origin/HEAD`, named for the change.
2. **Orient the agent** — the requirement, the section of the specification, the
   existing code it has to fit.
3. **Ask for a plan and the tests before any implementation**, and read the plan. A
   wrong approach costs a sentence to correct here and a rewrite to correct later.
4. **Implement in small steps**, running the tests after each one.
5. **Read every line of the diff.** A change you cannot explain is not reviewed, and it
   does not matter that it works.
6. **Commit with a message that says why**, and note in the pull request what the agent
   was asked to do.
7. **Open the pull request** once the tests pass locally.
8. **Close the session by updating the memory file** — see
   [Memory and session continuity](#memory-and-session-continuity).

### What the human owns

The agent instructions above constrain the agent. These are the reviewer's, and no
instruction given to an agent can substitute for them:

- **Review agent-written code as code from a stranger**, not as your own. 
- **Check the tests hardest.** An agent that wrote both the implementation and its
  test can agree with itself, and a test written after the fact asserts what the code
  does rather than what it should do.
- **Check every number and every citation** that will leave the repository. Fabricated
  references and numbers can look exactly like real ones.
- **Do not merge what you cannot explain.** If you could not defend the change in a
  group meeting, it is not ready. The fix is to understand it.

## Additional engineering conventions

- **Red-green TDD for bug fixes and regressions.** Write the failing test first; a fix
  with no test that failed before it is a guess.
- **CI enforces `clang-format`** on pull requests, against the `.clang-format` in the
  repository root. Format before pushing rather than in a follow-up commit.
- **Build and test** with `cmake -B build && cmake --build build && ctest --test-dir
  build`. The build needs C++17, MPI, Cabana with MPI and grid support, Silo, and
  GoogleTest; HYPRE is optional. Set `MPIEXEC_PREFLAGS` when MPI needs extra
  arguments, as `configs/github/host-config.cmake` does for OpenMPI.
- **Run every test at several rank counts**, and encode global, not local, indices in
  any value that crosses ranks. Tests that passed only on one rank hid several of the
  bugs `FUTURE.md` records.
- **Keep the halo deep enough for every stencil.** `ProblemManager.hpp` explains the
  current depth; a stencil that reads past it reads garbage silently, so check new
  kernels under a Kokkos debug build with bounds checking.
- **HYPRE built with CUDA breaks the serial backend**, so CUDA builds must run with
  `-p cuda`.
- **Determinism.** Same inputs, same outputs. Seed anything sampled and record the
  seed.
- **Inline single-use variables and functions**, prefer self-documenting code to
  explanatory comments, and follow YAGNI and DRY.

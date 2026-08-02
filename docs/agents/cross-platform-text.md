# Cross-platform text and fixture identity

All tracked text in this repository uses LF line endings. The root
`.gitattributes` enforces this at checkout on Linux, macOS, and Windows while
leaving files detected as binary untouched.

## Why this is a repository contract

Oracle fixtures, golden files, snapshots, and other independently generated
artifacts may be verified by a byte-level hash. Line endings are bytes: a
Windows checkout that silently changes LF to CRLF changes the artifact identity
even when its parsed contents are identical. This has repeatedly caused
Windows-only CI failures after Linux and macOS passed.

The durable fix belongs at the Git checkout boundary. Do not work around a
failure by accepting platform-specific hashes or normalizing the artifact inside
its identity verifier; either approach weakens the byte-identity contract.

## Adding byte-identified text artifacts

1. Store the artifact as UTF-8 with LF line endings.
2. Confirm the checkout rule with `git check-attr text eol -- <path>`; text
   artifacts must report `text: auto` and `eol: lf`.
3. Compute and record the hash from the staged Git content, not an editor's
   unstaged working copy: `git show :<path> | sha256sum`.
4. Keep at least one Windows CI target whenever a test uses `include_str!`, a
   snapshot, or a byte-level fixture hash.

If a file format genuinely requires different line endings, add the narrowest
possible override to `.gitattributes` and document why beside that rule.

# Shared code

Every code file of this pipeline is an unchanged copy of the file of the
same name in v_dynamic, the lab's dynamic pipeline (`../v_dynamic`),
whose static mode (`S.dynamic = false`) is this pipeline. A change to a
shared file is made in v_dynamic, tested there, and copied here:

```bash
bash tools/sync_from_v_dynamic.sh ../v_dynamic
```

It copies every shared file and rewrites `tools/v_dynamic_files.txt` (the v_dynamic
commit and each file's SHA-256). `tests/test_same_code_as_v_dynamic.m` fails
when a shared file was edited here or differs from v_dynamic's, when a code file
is neither shared nor this folder's own, or when the code here calls a
function that only v_dynamic has.

This folder's own files: `pipeline_settings.m`, `examples/`,
`tests/test_template.m`, `tests/test_same_code_as_v_dynamic.m`,
`tools/sync_from_v_dynamic.sh`, `README.md`, `docs/methods.md`,
`docs/validation.md`, this file, `CHANGES.md`, `LICENSE`, `LICENSE.md`, `.github/`,
`.gitattributes`, `.gitignore`.

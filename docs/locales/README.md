# Translation catalogs

This directory stores Sphinx i18n catalogs for PySME documentation.

Recommended workflow:

```bash
sphinx-build -b gettext docs docs/_build/gettext
sphinx-intl update -p docs/_build/gettext -l zh_CN
sphinx-build -D language=zh_CN -b html docs docs/_build/html-zh_CN
```

English files under `docs/` remain the single source of truth.
Translated content should be maintained under `docs/locales/<language>/LC_MESSAGES/`.

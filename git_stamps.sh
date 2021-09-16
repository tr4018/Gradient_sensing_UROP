#!/bin/sh

HEAD=$(git rev-parse HEAD)


printf 'static char *git_version_string=\"# Info: HEAD='${HEAD}'\\n\"\n'
# https://stackoverflow.com/questions/2179722/checking-out-old-file-with-original-create-modified-timestamps
for FILE in $(git ls-files)
do
    TIME=$(git log --pretty=format:%cd -n 1 --date=iso -- "$FILE")
    printf '\"# Info: %s :: %s\\n\"\n' "$FILE" "$TIME"
done
echo ';'

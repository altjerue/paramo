#!/bin/bash
echo -e "Commit?"
read com
git add -A
git status
git commit -m "$com"
git push https://zkd:wmkwwk11@github.rcac.purdue.edu/JesusRuedaBecerrilGroup/Paramo.git --all

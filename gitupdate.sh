#!/bin/bash
echo -e "Commit?"
read com
git add -A
git status
git commit -m "$com"
git push origin https://zkd:wmkwwk11@github.rcac.purdue.edu/JesusReudaBecerrilGroup/Paramo.git

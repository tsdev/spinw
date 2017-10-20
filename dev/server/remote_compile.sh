# push changes
cd ~/spinw_git
git add .
git commit -a -m "Minor bug fix for remote testing"
git push
# login in a separate window before to make passwordless login
ssh nemesis './compile.sh'

Last login: Thu Jun  2 20:27:27 on ttys002
cd%                                                                                                     simondriscoll@reader172-210 ~ % cd icepack-dirs/runs/mycase 
simondriscoll@reader172-210 mycase % ls
compile				ice_diag.icefree		icepack.runlog.220602-201146
env.conda_macos			ice_diag.land			icepack.settings
history				ice_diag.slab			icepack_in
ice_diag.full_ITD		icepack				restart
simondriscoll@reader172-210 mycase % rm -fr *
zsh: sure you want to delete all 12 files in /Users/simondriscoll/icepack-dirs/runs/mycase [yn]? y
simondriscoll@reader172-210 mycase % ls
simondriscoll@reader172-210 mycase % source ~/.profile 
simondriscoll@reader172-210 mycase % ls -rlt
total 0
simondriscoll@reader172-210 mycase % pwd
/Users/simondriscoll/icepack-dirs/runs/mycase
simondriscoll@reader172-210 mycase % cd
simondriscoll@reader172-210 ~ % cd Icepack
cd: no such file or directory: Icepack
simondriscoll@reader172-210 ~ % 
simondriscoll@reader172-210 ~ % 
simondriscoll@reader172-210 ~ % 
simondriscoll@reader172-210 ~ % 
simondriscoll@reader172-210 ~ % 
simondriscoll@reader172-210 ~ % git clone https://github.com/simondriscoll/Icepack.git
Cloning into 'Icepack'...
remote: Enumerating objects: 5119, done.
remote: Counting objects: 100% (318/318), done.
remote: Compressing objects: 100% (185/185), done.
remote: Total 5119 (delta 177), reused 238 (delta 128), pack-reused 4801
Receiving objects: 100% (5119/5119), 8.49 MiB | 6.17 MiB/s, done.
Resolving deltas: 100% (3614/3614), done.
simondriscoll@reader172-210 ~ % git status
fatal: not a git repository (or any of the parent directories): .git
simondriscoll@reader172-210 ~ % cd Icepack 
simondriscoll@reader172-210 Icepack % git status
On branch main
Your branch is up to date with 'origin/main'.

nothing to commit, working tree clean
simondriscoll@reader172-210 Icepack % 
simondriscoll@reader172-210 Icepack % 
simondriscoll@reader172-210 Icepack % 
simondriscoll@reader172-210 Icepack % 
simondriscoll@reader172-210 Icepack % git status
On branch main
Your branch is up to date with 'origin/main'.

Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git restore <file>..." to discard changes in working directory)
	modified:   columnphysics/icepack_therm_vertical.F90
	modified:   configuration/driver/icedrv_diagnostics.F90
	modified:   configuration/scripts/machines/env.conda_macos

no changes added to commit (use "git add" and/or "git commit -a")
simondriscoll@reader172-210 Icepack % git add -A
simondriscoll@reader172-210 Icepack % git commit -m "Added the standard changes to columnphysics/icepack_therm_vertical.F90, configuration/driver/icedrv_diagnostics.F90, configuration/scripts/machines/env.conda_macos."
[main 6519d92] Added the standard changes to columnphysics/icepack_therm_vertical.F90, configuration/driver/icedrv_diagnostics.F90, configuration/scripts/machines/env.conda_macos.
 Committer: SimonDriscoll <simondriscoll@reader172-210.ouls-open.ox.ac.uk>
Your name and email address were configured automatically based
on your username and hostname. Please check that they are accurate.
You can suppress this message by setting them explicitly. Run the
following command and follow the instructions in your editor to edit
your configuration file:

    git config --global --edit

After doing this, you may fix the identity used for this commit with:

    git commit --amend --reset-author

 3 files changed, 113 insertions(+), 2 deletions(-)
simondriscoll@reader172-210 Icepack % 


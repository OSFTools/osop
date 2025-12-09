# Initial set up

### 2.  **Basic set up:**

For common users of code or developers within code frameworks, this
guide may cover redundant ideals. As such a brief overview of what is
included within this section is as follows; download and set up for a
package management system (Conda, in this case supported by
Miniforge), download and set up of Gitbash, mirroring of a repository
and finally environment installation. It would be expected that most
common coders have some or all of the considerations covered within
this section and as such it would be advised to move ahead. If issues
later present in installation and running of the OSOP toolkit it may
be pertinent to recheck this section for guidance.

### **2a. Installation of Miniforge:**

Conda is an open-source packet management system that is used to
install, maintain and handle software tools generally falling within
the realm of coding spaces. Miniforge is a lightweight installer for
Conda that is community-lead and openly licensed. It is not
necessarily the most user-friendly software; however, its minimalist
system ensures that it doesn't requisite heavily on system
requirements and avoids clutter from larger installers like Anaconda.
To be able to set up and get full functionality of the OSOP toolkit, a
package management software is needed to set up the environments
needed for the toolkit to run as well as maintain the modules long
term.

Guidance for setting up Miniforge can be found here:
[conda-forge/miniforge: A conda-forge
distribution.](https://github.com/conda-forge/miniforge)

It would be recommended here to follow through with the advice laid
out by the Miniforge developers for your specific situation. Miniforge
is the underpinning software that allows functionality of everything
else to come, so taking time here to get things smoothly working is
advised.

### **2b. Installation of Git:**

Git is a version control system that allows for multiple developers to
work on one project simultaneously without interference. Meanwhile,
Gitbash is emulator system that allows for communication with Git
services whilst enabling functionality of Unix based features. This
set up is based in two parts. One is the installation, and the other
is integration with Miniforge.

### **2b.i. Installation:**

It is important that during the installation of Git, Gitbash is
selected within this process. Otherwise set up systems are up to the
user's discretion.

Download of Git can be found here: [Git -
Downloads](https://git-scm.com/downloads)

Guidance for setting up Git can be found here: [Git Guides - install
git](https://github.com/git-guides/install-git)

Git is the main system for later downloading and mirroring the
repository that contains the OSOP toolkit. It is also the "missing
piece" in allowing for full functionality of the product, issues with
running the baseline product -- especially on windows system, may be
traced back here.

### **2b.ii. Cross functionality with Miniforge:**

Git does not necessarily function with Conda out the gate. Before
setting up the toolkit Git needs to be configured with Conda so that
packages can be loaded within the Git space.

The first check here is to see if its configured to work from
download. Open a Gitbash terminal via the search bar at the bottom --
type Git bash and it should come up. Upon completion a terminal should
appear. Typing `conda --version` will result in two options. The first
means that Git is configured with conda straight away, this will be if
it displays something along the lines of conda 25.3.1. The other
option is that it returns conda command not found, this means that
conda is not added to the systems PATH.

To fix this issue the following is advised. Type `nano ~/.bashrc` and
enter onwards. A text file system should appear. At the bottom of this
document (it may be blank upon opening, just type from there) place in
the following lines:

```
. "/c/Users/\<username\>/miniforge3/etc/profile.d/conda.sh"

conda activate base
```

**note replace \<username\> with the user of the OS.**

From here press ctrl + x, followed by y and then the enter key. Exit
Gitbash, open it again and type `conda --version`, this should now
work and return a conda 25.3.1 or similar.

This is the final step in software downloads that will underpin the
OSOP toolkit. Most of what is outlined here is necessary to run
python-based repositories from GitHub. Here onwards is about settling
in the OSOP product before finally being able to customise it to your
specific usage case.

### **2c. Mirroring the repository or Downloading:**

The baseline guidance here is to download the repository as a zip file
and extract. To do this go to this link [OSFTools/osop: Objective
Seasonal Outlooks Package](https://github.com/OSFTools/osop/) and
press on the green button that says code. Follow through with download
as ZIP on the drop-down box that appears. Once downloaded open up
Gitbash, and navigate to the download. This will usually be a case of
typing cd Downloads + enter, then ls and enter to confirm that you can
see a osop-main.zip file. From here type `extract osop-main.zip` and
enter onwards. This should create a osop-main/ file that can be looked
into and no longer has the .zip suffix.

Alternatively you can mirror the repository. Mirroring a repository is
what is generally considered to be a shallow copy of the original.
This is to mean that it copies the system locally but without certain
abilities in place. Most base level users should find that mirroring
is plenty for straight forward level usage and manipulation of the
code. However, developers or users interested in furthering the code
may want to consider forking the repo instead.

Guidance for mirroring a repository can be found here: [Duplicating a
repository - GitHub
Docs](https://docs.github.com/en/repositories/creating-and-managing-repositories/duplicating-a-repository?r292)

The specific repo can be found here: [OSFTools/osop: Objective
Seasonal Outlooks Package](https://github.com/OSFTools/osop/)

At this stage the repository is now handleable locally. The user
should take a moment to gather bearings and explore the functionality
of what has been done so far. Check that Conda is accessible within
Gitbash, and that tour around the repository to get a gist of layout
and location.

### **2d. Installing the environment:**

The final step in running the OSOP toolkit is to set up the
dependencies of the python code within the package. The dependencies
are modules that allow python code to function. They are essentially
shorthand for larger code that has been written and saved out by a
developer to allow for more code to be written that utilises these
base tools. The dependency of a project is usually included in the
"environment.yml" file found in most repositories. The first step will
be navigating within the osop-main file discussed in section 2c. This
is usually a case of `cd osop-main` within the directory that holds
this. If following the download zip section then this will be in
downloads. Typing `ls` will bring up a list of the files contained
here. One of which should be the environment.yml

Guidance for installing and activating an environment from a yaml file
can be found here: [Managing environments --- conda 25.7.1.dev65
documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
under the "*Creating an environment from an environment.yml file*"
subheading.

This concludes the set-up process of the OSOP toolkit. The next
section is to outline running the toolkit and then finally exploring
how to get the best usage out of it for your specific needs.

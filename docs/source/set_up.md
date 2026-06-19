# Initial set up

This section covers basic concepts and so experienced Python users can
quickly move to the later parts. If issues are found later in
installation and running of the OSOP toolkit check this section for
guidance.

A brief overview of what is
included within this section is as follows: 

* download and set up a package management system (Conda, in this case supported by
Miniforge)
* download and set up of Gitbash (for Window users only)
* mirroring of a repository
* environment set up

## Installation of Miniforge

To be able to set up and get full functionality of the OSOP toolkit,
package management software is needed to set up the environments
needed for the toolkit to run as well as to maintain the modules long
term.

Conda is an open-source packet management system that is used to
install, maintain and handle software. Miniforge is a lightweight
installer for Conda that is community-lead and openly licensed. It is
primarily a command line based application (unlike Anaconda Navigator
which uses a Graphical User interface); however, its minimalist and
open system ensures that it can be used without licences even by
larger organisations and has lower system requirements than larger
installers like Anaconda.

Guidance for setting up Miniforge can be found here:
[conda-forge/miniforge: A conda-forge
distribution.](https://github.com/conda-forge/miniforge)

We recommend following through the advice of the Miniforge developers
for your specific situation. Miniforge is the underpinning software
that supports everything else, so taking time here to get things
working smoothly is advised.

## Installation of Git and Gitbash on Windows

This section is for Windows users.

[Git](https://git-scm.com/) is a version control system that allows for multiple developers to
work on one project simultaneously. [Gitbash](https://gitforwindows.org/) is
an emulator system for Windows 
that allows for communication with Git
services whilst enabling functionality of Unix based features. This
set up is based in two parts. One is the installation, and the other
is integration with Miniforge.

### Installation

It is important that during the installation of Git, Gitbash is
selected within this process. Otherwise set up systems are up to the
user's discretion.

Download of Git can be found here: [Git -
Downloads](https://git-scm.com/downloads)


Git is the main system for later downloading and mirroring the
repository that contains the OSOP toolkit. It is also the "missing
piece" in allowing for full functionality of the product, issues with
running the baseline product -- especially on windows system, may be
traced back here.

### Miniforge in Gitbash (Windows only)

Gitbash does not necessarily work with Conda straight away. Before
setting up the toolkit Gitbash needs to be configured so that conda
packages can be loaded within the Git space.

The first check here is to see if its configured to work from
download. Open a Gitbash terminal via the search bar at the bottom --
type Git bash and it should come up. Upon completion a terminal should
appear. Typing `conda --version` will result in two options. If this
displays something like `conda 25.3.1` thenthis  means that it is
configured with conda and you can go onto the next section. The other
option is that it returns `bash: conda: command not found`, this means
that conda is not added to the systems PATH.

To fix this issue the following is advised. Type `nano ~/.bashrc` and
press enter. A text file system should appear. At the bottom of this
document (it may be blank upon opening, just type from there) place in
the following lines:

```bash
. "/c/Users/\<username\>/miniforge3/etc/profile.d/conda.sh"
conda activate base
```

**note replace \<username\> with the user of the OS.**

From here press ctrl + x, followed by y and then the enter key. Exit
Gitbash, open it again and type `conda --version`, this should now
work and return  `conda 25.3.1` or similar.

This is the final step in software downloads that will underpin the
OSOP toolkit. Most of what is outlined here is necessary to run
python-based repositories from GitHub. Here onwards is about settling
in the OSOP product before finally being able to customise it to your
specific usage case.

## Install Git (Linux)
Guidance for setting up Git on Linux can be found here: [Git Guides -
install git](https://github.com/git-guides/install-git). 

## Getting the OSOP code

The easiest way to get start with the code is to download the repository as a zip file
and extract. To do this go to this link [OSFTools/osop: Objective
Seasonal Outlooks Package](https://github.com/OSFTools/osop/) and
press on the green button that says code. Follow through with download
as ZIP on the drop-down box that appears. Once downloaded open up
Gitbash, and navigate to the download by
typing `cd Downloads` + enter, then `ls` and enter to confirm that you can
see a osop-main.zip file. From here type `unzip osop-main.zip` and
enter. This should create a osop-main/ folder with the code in it.

To move this to your home directory type:
`mv osop-main $HOME`

At this stage the repository is now usable locally. The user
should take a moment to browse the files in the code and explore the functionality of the code.

## Mirror the code
Alternatively you can mirror the repository. Mirroring a repository is
what is generally considered to be a shallow copy of the original.
This means that it copies the system locally but without certain
abilities. Most base level users should find that mirroring
is good enough for straightforward usage and manipulation of the
code. However, developers or users interested in improving the code
may want to consider cloning the repo instead.

Guidance for mirroring a repository can be found here: [Duplicating a
repository - GitHub
Docs](https://docs.github.com/en/repositories/creating-and-managing-repositories/duplicating-a-repository?r292)

The specific repo can be found here: [OSFTools/osop: Objective
Seasonal Outlooks Package](https://github.com/OSFTools/osop/)

## Installing the environment

The final step in running the OSOP toolkit is to set up the Python 
code to run the package. This makes available the modules that allow the Python code to function. 

The first step will be navigating to the osop-main folder discussed in the
previous section. This
is usually a case of `cd osop-main` within the directory that holds
this. If following the download zip section then this will be in
downloads. Typing `ls` will bring up a list of the files contained
here. One of which should be the file `osop-spec-file.yml`.
This contains a detailed list of the packages and the versions 
needed to run the code.

Guidance for installing and activating an environment from a yaml file
can be found here: [Managing environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)
under the "*Creating an environment from an environment.yml file*"
subheading.

This concludes the set-up process of the OSOP toolkit. The next
section outlines running the toolkit and finally exploring
how to get the best use out of it for your specific needs.

## Create CDS account and set up for seasonal downloads

To be able to download data you will need to first sign up for
an ECMWF account if you do not have one already. Start by creating
an account from 
[Climate Data Store](https://cds.climate.copernicus.eu/#!/home). 
Click "Login-Register" on the top right and "Register New User".

Then follow the instructions on [CDSAPI setup - Climate Data
Store](https://cds.climate.copernicus.eu/how-to-api). This 
shows you how to get your CDS API key (like a password).

You will also need to visit the download page for ERA5 and Seasonal
datasets when logged in to the website and accept the terms 
and conditions of use.
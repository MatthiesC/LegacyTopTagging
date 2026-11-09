PROMPT="%{$fg_bold[blue]%}%* %{$fg_bold[green]%}[%n@%m]%{$fg[blue]%}%0~%# %{$reset_color%}"
#PS1='\[\033[1;34;40m[\033[1;31;40m\u@\h:\w\033[1;34;40m]\033[1;37;40m $\033[0;37;0m\] ' # old prompt used since beginning of master thesis until Nov 1, 2019
PS1='[ \[\033[1;34m\]\u\[\033[0m\]@\[\033[1;33m\]\h\[\033[0m\]:\[\033[1;34m\]\w\[\033[0m\] ]\n $ '
#➤ ➤ ➤  '

alias ls='ls --color=auto'
alias la='ls -a'
alias lh='ls -lh'
alias ll='ls -lh'

alias proxy='voms-proxy-init --voms cms:/cms/dcms'

# need to be on the NAF to use these:
alias print='lpr -H uh2usitcups -P 68_130_a4_bw_km'
alias print_color='lpr -H uh2usitcups -P 68_130_a4_cl_km'

# for access to Pixel DCS:
alias windoostunnel="ssh -Y -l cmatthie -C -L 9999:cerntscms.cern.ch:3389 lxplus.cern.ch"
alias windoosdesk="rdesktop -g 90% -P -z -x l -u cmatthie -d CERN localhost:9999"

# set up Anaconda2 environment:
alias ana2='
export PATH=/data/dust/user/matthies/anaconda2/bin:$PATH
source activate /data/dust/user/matthies/anaconda2/envs/myenv
'

alias dusthome='cd /data/dust/user/matthies/'

export PATH="/data/dust/user/matthies/texlive/2022/bin/x86_64-linux/:$PATH"

alias naf11='ssh naf-cms11'
alias naf12='ssh naf-cms12'
alias naf13='ssh naf-cms13'
alias naf14='ssh naf-cms14'
alias naf15='ssh naf-cms15'
alias naf16='ssh naf-cms16'
alias naf17='ssh naf-cms17'
alias naf18='ssh naf-cms18'

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/dust/user/matthies/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/dust/user/matthies/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/data/dust/user/matthies/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data/dust/user/matthies/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

#conda activate ana3
alias pyconda3='/data/dust/user/matthies/anaconda3/envs/ana3/bin/python'
alias ipyconda3='/data/dust/user/matthies/anaconda3/envs/ana3/bin/ipython'

alias sourcecms="source /cvmfs/cms.cern.ch/cmsset_default.sh"
#alias el7="cmssw-el7 --contain --bind /afs:/afs --bind /nfs:/nfs --bind /pnfs:/pnfs --bind /cvmfs:/cvmfs --bind /var/lib/condor:/var/lib/condor"
alias el7="cmssw-el7 --contain --bind /afs:/afs --bind /data:/data --bind /pnfs:/pnfs --bind /cvmfs:/cvmfs --bind /var/lib/condor:/var/lib/condor"
alias sel7="sourcecms; el7"

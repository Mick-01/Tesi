FROM opencfd/openfoam-default:latest

# Install nano, paraview, and some extras for usability
RUN apt-get update && apt-get install -y --no-install-recommends \
      nano \
      bash-completion \
      coreutils \
    && rm -rf /var/lib/apt/lists/*

# Set a colored bash prompt (username@host in green, path in blue, $ in yellow)
RUN echo '' >> /etc/bash.bashrc && \
    echo '# Custom colored prompt' >> /etc/bash.bashrc && \
    echo 'PS1="\[\033[01;32m\]\u@\h \[\033[01;34m\]\w \[\033[01;33m\]\$\[\033[00m\] "' >> /etc/bash.bashrc && \
    echo 'alias ls="ls --color=auto"' >> /etc/bash.bashrc && \
    echo 'alias grep="grep --color=auto"' >> /etc/bash.bashrc


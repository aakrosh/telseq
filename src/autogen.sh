#!/bin/sh

set -ex
aclocal
autoconf -Wall
autoheader
automake -a


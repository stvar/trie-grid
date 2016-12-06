#!/bin/bash

# Copyright (C) 2016  Stefan Vargyas
# 
# This file is part of Trie-Grid.
# 
# Trie-Grid is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Trie-Grid is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Trie-Grid.  If not, see <http://www.gnu.org/licenses/>.

#
# File generated by a command like:
# $ gen-test -T unique
#

[[ "$1" =~ ^-u[0-9]+$ ]] &&
u="${1:2}" ||
u=""

diff -u$u -L unique.old <(echo \
'$ test -x ./trie
$ print() { printf '\''%s\n'\'' "$@"; }
$ trie() { ./trie -C -u >/dev/null; }
$ print b a|trie
trie: error: line #2: input is not sorted: "a"
command failed: print b a|trie
$ print a a b|trie
trie: error: line #2: input is not uniquely sorted: "a"
command failed: print a a b|trie
$ print a b a|trie
trie: error: line #3: input is not sorted: "a"
command failed: print a b a|trie
$ print a b b|trie
trie: error: line #3: input is not uniquely sorted: "b"
command failed: print a b b|trie
$ print b b a|trie
trie: error: line #3: input is not sorted: "a"
command failed: print b b a|trie
$ print a A|trie
trie: error: line #2: input is not sorted: "A"
command failed: print a A|trie'
) -L unique.new <(
echo '$ test -x ./trie'
test -x ./trie 2>&1 ||
echo 'command failed: test -x ./trie'

echo '$ print() { printf '\''%s\n'\'' "$@"; }'
print() { printf '%s\n' "$@"; } 2>&1 ||
echo 'command failed: print() { printf '\''%s\n'\'' "$@"; }'

echo '$ trie() { ./trie -C -u >/dev/null; }'
trie() { ./trie -C -u >/dev/null; } 2>&1 ||
echo 'command failed: trie() { ./trie -C -u >/dev/null; }'

echo '$ print b a|trie'
print b a|trie 2>&1 ||
echo 'command failed: print b a|trie'

echo '$ print a a b|trie'
print a a b|trie 2>&1 ||
echo 'command failed: print a a b|trie'

echo '$ print a b a|trie'
print a b a|trie 2>&1 ||
echo 'command failed: print a b a|trie'

echo '$ print a b b|trie'
print a b b|trie 2>&1 ||
echo 'command failed: print a b b|trie'

echo '$ print b b a|trie'
print b b a|trie 2>&1 ||
echo 'command failed: print b b a|trie'

echo '$ print a A|trie'
print a A|trie 2>&1 ||
echo 'command failed: print a A|trie'
)


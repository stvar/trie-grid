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
# $ gen-test -T base
#

[[ "$1" =~ ^-u[0-9]+$ ]] &&
u="${1:2}" ||
u=""

diff -u$u -L base.old <(echo \
'$ test -x ./trie
$ print() { printf '\''%s\n'\'' "$@"; }
$ print aabc aabd xxx xyz a|sort -u|./trie -u -c
{
    "a": "a" {
        "ab" {
            "c": "aabc"
            "d": "aabd"
        }
    }
    "x" {
        "xx": "xxx"
        "yz": "xyz"
    }
}
$ print aabc aabd xxx xyz a|shuf|./trie -U -c
{
    "a": "a" {
        "ab" {
            "c": "aabc"
            "d": "aabd"
        }
    }
    "x" {
        "xx": "xxx"
        "yz": "xyz"
    }
}
$ print aabc aabd xxx xyz a|sort -u|./trie -u -w
{
    '\''a'\'': "a" {
        '\''a'\'' {
            '\''b'\'' {
                '\''c'\'': "aabc"
                '\''d'\'': "aabd"
            }
        }
    }
    '\''x'\'' {
        '\''x'\'' {
            '\''x'\'': "xxx"
        }
        '\''y'\'' {
            '\''z'\'': "xyz"
        }
    }
}
$ print aabc aabd xxx xyz a|shuf|./trie -U -w
{
    '\''a'\'': "a" {
        '\''a'\'' {
            '\''b'\'' {
                '\''c'\'': "aabc"
                '\''d'\'': "aabd"
            }
        }
    }
    '\''x'\'' {
        '\''x'\'' {
            '\''x'\'': "xxx"
        }
        '\''y'\'' {
            '\''z'\'': "xyz"
        }
    }
}'
) -L base.new <(
echo '$ test -x ./trie'
test -x ./trie 2>&1 ||
echo 'command failed: test -x ./trie'

echo '$ print() { printf '\''%s\n'\'' "$@"; }'
print() { printf '%s\n' "$@"; } 2>&1 ||
echo 'command failed: print() { printf '\''%s\n'\'' "$@"; }'

echo '$ print aabc aabd xxx xyz a|sort -u|./trie -u -c'
print aabc aabd xxx xyz a|sort -u|./trie -u -c 2>&1 ||
echo 'command failed: print aabc aabd xxx xyz a|sort -u|./trie -u -c'

echo '$ print aabc aabd xxx xyz a|shuf|./trie -U -c'
print aabc aabd xxx xyz a|shuf|./trie -U -c 2>&1 ||
echo 'command failed: print aabc aabd xxx xyz a|shuf|./trie -U -c'

echo '$ print aabc aabd xxx xyz a|sort -u|./trie -u -w'
print aabc aabd xxx xyz a|sort -u|./trie -u -w 2>&1 ||
echo 'command failed: print aabc aabd xxx xyz a|sort -u|./trie -u -w'

echo '$ print aabc aabd xxx xyz a|shuf|./trie -U -w'
print aabc aabd xxx xyz a|shuf|./trie -U -w 2>&1 ||
echo 'command failed: print aabc aabd xxx xyz a|shuf|./trie -U -w'
)


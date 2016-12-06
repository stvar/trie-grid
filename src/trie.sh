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

trie-gcc-version()
{
    local p="${1:-gcc}"

    set -o pipefail
    "$p" --version|
    sed -nr '
        1 {
            s/^gcc.*\s([0-9]+\.[0-9]+\.[0-9]+)\s.*$/\1/p
            q
        }'|
    awk -F. '{
        printf("%d%02d%02d\n", $1, $2, $3)
    }'
}


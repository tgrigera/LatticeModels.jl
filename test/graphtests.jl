# graphtests.jl -- tests for graphs and lattices
#
# Copyright (C) 2023 Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For details, see the file LICENSE in the root directory, or
# check <https://www.gnu.org/licenses/>.

function SQLattice_test()
    @testset "Square lattice tests" begin
        L=5
        M=7
        lat=SQLattice_open{Tuple{Int,Int}}(L,M)
        for j=1:M,i=1:L
            lat[i,j]=(i,j)
        end
        n=Ref(0)
        foreach_bond(lat) do i,j
            n[]+=1
        end
        @test n[] == 2*L*M-L-M

        rs=random_site(lat)
        foreach_neighbour!(rs) do s
            s[]=(0,0)
        end
        i,j=rs.I[1],rs.I[2]
        @test (i>=L || lat[i+1,j]==(0,0)) &&
            (i<=1 || lat[i-1,j]==(0,0)) &&
            (j>=M || lat[i,j+1]==(0,0)) &&
            (j<=1 || lat[i,j-1]==(0,0))

        lat=SQLattice_periodic{Tuple{Int,Int}}(L,M)
        for j=1:M,i=1:L
            lat[i,j]=(i,j)
        end
        n=Ref(0)
        foreach_bond(lat) do i,j
            n[]+=1
        end
        @test n[] == 2*L*M
        rs=random_site(lat)
        foreach_neighbour!(rs) do s
            s[]=(0,0)
        end
        i,j=rs.I[1],rs.I[2]
        @test (i>=L || lat[i+1,j]==(0,0)) &&
            (i<L || lat[1,j]==(0,0)) &&
            (i<=1 || lat[i-1,j]==(0,0)) &&
            (i>1 || lat[L,j]==(0,0)) &&
            (j>=M || lat[i,j+1]==(0,0)) &&
            (j<M || lat[i,1]==(0,0)) &&
            (j<=1 || lat[i,j-1]==(0,0)) &&
            (j>1 || lat[i,M]==(0,0))

    end

end

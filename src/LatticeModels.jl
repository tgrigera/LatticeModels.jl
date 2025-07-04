# LatticeModels.jl
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

module LatticeModels

include("./graphs/Graphs.jl")

include("./aux/LatticeTrajectory.jl")

include("./models/IsingModel.jl")

include("./models/Sandpiles.jl")

include("./models/Heisenberg.jl")

end

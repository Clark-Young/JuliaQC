# Description: This Julia Code is used to interpret the input file for quantum chemistry calculation.
# We accept the input file in the format of .gjf file (Gaussian input file) or .xyz file (xyz format).
# If you use .xyz file, you will be prompted to input the charge and multiplicity of the molecule (if not specified in the second line), the method and basis set you want to use.
# You need know that the expansion must match the format of the input file.
# Author: Chunyu Yang
# Email: chunyu.yang@duke.edu, chyyangustc@outlook.com (permanent)
# Organization: Department of Chemistry, Duke University
# Date: 2024-12-13
# Reference: 


module InputIntrepreter

struct illegalInputError <: Exception 
    msg::String
end

function callIntrepreter(filePath::String)
    dotindex = findlast(==('.'), filePath)
    if isnothing(dotindex)
        throw(illegalInputError("No extension found."))
    else
        extension = filePath[dotindex+1:end]
        if extension == "gjf"
            gjfIntrepreter(filePath)
        elseif extension == "xyz"
            xyzIntrepreter(filePath)
        else
            throw(illegalInputError("The extension ."+extension+" is not supported. Please use .gjf or .xyz file."))
        end
    end
end

function xyzIntrepreter(filePath::String)
    # This function is used to interpret the xyz file.
    # The xyz file is a simple format that contains the number of atoms, and the coordinates of each atom.
    # The first line is the number of atoms, and the second line is the comment, here it could be the charge and multiplicity of the molecule.
    # The following lines are the coordinates of each atom.
    # The format is as follows:
    # 3
    # 0 1 H2O # The first two number are the charge and multiplicity of the molecule.
    # O 0.0 0.0 0.0
    # H 0.0 0.0 1.0
    # H 0.0 1.0 0.0
    # The function will return the atomlist, charge, multiplicity, method, and basis set.

    
end


function gjfIntrepreter(filePath::String)
    # This function is used to interpret the Gaussian input file.
    # You can learn about the Gaussian input file from the Gaussian manual. If you know little about it, use the .xyz file.


    
end

end


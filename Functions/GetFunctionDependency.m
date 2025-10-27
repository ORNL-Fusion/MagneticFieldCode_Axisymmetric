function [fList,pList] = GetFunctionDependency(scriptName)
% GetFunctionDependency
% -------------------------------------------------------------------------
% This function takes the name of a matlab script/function and determines
% which other function are needed to succesfull run the code.
% This provides a list of the "children" functions.
% This is useful in determining which functions need to be fetched to make
% a script work correctly in case it is shared with another person
% -------------------------------------------------------------------------
% INPUT:
% -------------------------------------------------------------------------
% scriptName: string that contains the name of the script we want to
% determine which other functions and scripts are need to run it.
% -------------------------------------------------------------------------
% OUTPUT:
% fList: cell that contains the paths to each of the functions or scripts
% required to run (scriptName).m. Each cell contains a string address.
% pList: Returns information about the toolboxes required to run the script
% -------------------------------------------------------------------------

[fList,pList] = matlab.codetools.requiredFilesAndProducts(scriptName);

end


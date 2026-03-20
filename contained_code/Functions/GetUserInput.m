function [answerDouble] = GetUserInput(inputStructure)
% Fields that can be defined by the user:
% prompt
% name
% defaultAnswer

% Prompt
prompt = inputStructure.prompt; 

% Name
try
    name = inputStructure.name;
catch
    name = '';
end

% Default answer
try 
    defaultAnswer = inputStructure.defaultAnswer;
catch
    defaultAnswer = {'1'};
end

% Window style
options.WindowStyle = 'normal';

% Get answer from user
try % prompt is a cell
    answer = inputdlg(prompt,name,[1,length(prompt{:})+30],defaultAnswer,options);
catch % prompt is an array
    answer = inputdlg(prompt,name,[1,length(prompt)+30],defaultAnswer,options);
end

% Interprete input data
if isempty(answer); 
    answer{1} = '0';
elseif strcmp(answer{:},'3')
    % kill program
end

% convert answer into a double type
answerDouble = str2num(answer{:});

end


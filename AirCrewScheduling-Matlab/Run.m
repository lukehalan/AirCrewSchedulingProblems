sppnw = "sppnw43.txt";
maximumIteration = 1000;
Start(sppnw, maximumIteration)

% If you want test 30 Independetn run please uncomment below codes
% and comment Start function from above

% results = [];
% repetition = 30;
% for i = 1:repetition
% disp(['Run #: ' , num2str(i)]);    
% [C , S, V] = Start(sppnw, maximumIteration);
% results = [results; [C,V]];
% end
% 
% disp(char(10));
% disp(['Results after ' , num2str(repetition) , ' times run']);
% disp(results);
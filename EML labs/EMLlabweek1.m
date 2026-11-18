
Page
1
of 5
printf("============================================================\n")
printf("Intro to Octave\n")
printf("Kyle Beggs\n")
display(date())
printf("============================================================\n\n")
diary output.txt  % save the output to a file with this name
diary on
% this is a comment
printf("This is how to print to the screen.\n")
% Don't worry about learning all of this right away. We will be covering the
% concepts you need to know within each lab to complete that lab.
printf("\n\nDefining variables\n")
printf("==================\n\n")
% all variables default to double types even though it prints as an int.
x = 2
% add a semi-colon ; to the end to suppress output
x = 3;
% define an array
x = [1, 2, 3]
% or you can leave out the commas ,
%  (but I like the commas as it makes it easier to read)
x = [1 2 3];
% Notice both are sideways - they are ROW vectors. Octave distinguishes
%  from ROW vectors and COLUMN vectors as we do in formal mathematics.
% To write a COLUMN vector, use semi-colons ;
x = [4; 5; 6]
% To create a matrix we combine these ideas of creating row and column vectors.
A = [1,2,3; 4,5,6; 7,8,9]
printf("\n\nAccessing elements and slicing\n")
printf("==============================\n\n")
% to access a particular element in the array/vector
x(2)
% if we don't assign it to a variable, Octave auto assigns it to a variable
%  called 'ans' - notice this in the command window
% notice array indexing start at 1, not 0 like in C.
% extract an element in the matrix as
A(2,3)
% shorthand way to create an array with fixed steps
x = 1:5  % defaults to steps/increments of 1
% can specify the step size as the middle number
x = 0:0.25:1
% can extract a whole row or column, or even partial rows/columns
% extract the first column. : means take all elements in that direction
A(:,2)
% since first index refers to the row and second index is column, we are saying
% take all the rows of the second column
% we can also only take some of the rows
A(1:2,2)
% or
A(2:end,2)  % end keyword corresponds to the last element
printf("\n\nOperations\n")
printf("==========\n\n")
% add, substract, multiply, divide is obvious
a = 1;
b = 2;
z = a+b
z = a-b
z = a*b
z = a/b
% dot product
x = [1; 2; 3]  % column vector
y = [4, 5, 6]  % row vector
w = y*x
% w = x*y
% The above won't work because in formal mathemtics the dot product is 
%  a row vector into a column vector, not other way around.
% We can premultiply a matrix and vector to get another vector
b = A*x
% b = A*y won't work because wrong dimensions, same reason as above
printf("\n\nMemory\n")
printf("======\n\n")
% delete memory
clear x  % 'clear all' to delete all memory
% dynamically adding elements to an array
x = [1, 2]
x(3) = 3;
x
% notice the array is one element longer now because we set a number
%  to an index that was one more than the end of the array
% but this is NOT efficient. If we know how many elements are in the array
% we should allocate memory all at once using the zeros(rows,columns) function
x = zeros(1,3)
% then we will insert
x(1) = 1;
x(2) = 2;
x(3) = 3;
x
printf("\n\nLoops and if/else statements\n")
printf("============================\n\n")
% this should be familiar, just showing the syntax
% get how many elements are in the array
rows = size (A,1)  % the second arguemnt in size function tells which direction
columns = size (A,2)
% can use length for an array, doesn;t matter if it is a row or column vector
num_elements = length(x)
x = 1:5
for i = 1:length(x)
  x(i)
end
printf("\n")
i = 1;
while i <= 5
  x(i)
  ++i;
end
printf("\n")
x = 2;
y = 2;
z = 1;
if x==y
  printf("yeah buddy\n")
elseif x==z
  printf("weird.\n")
else
  printf("oh no\n")
end
printf("\n\nFunctions\n")
printf("=========\n\n")
% declare a function of a single input and output
function y = yourgrade(x)
  y = x+10;
end
yy = yourgrade(90)
printf("\n")
% define a function with multiple inputs and outputs
function [percent, letter] = wHaTisMygRAdE(quizzes,projects,exam1,exam2,final)
  percent = 0.2*mean(quizzes) + 0.3*mean(projects) + ...
    0.15*exam1 + 0.15*exam2 + 0.2*final;
  if percent >= 90
    letter = "A";
  elseif percent >= 88
    letter = "A-";
  elseif percent >= 86
    letter = "B+";
  elseif percent >= 80
    letter = "B";
  elseif percent >= 78
    letter = "B-";
  elseif percent >= 76
    letter = "C+";
  elseif percent >= 70
    letter = "C";
  else
    letter = "RIP";
  end
end
% another one for fun where you dont return anything
function what_do_i_need_on_my_final(target, quizzes, projects, exam1, exam2)
  current = 0.2*mean(quizzes) + 0.3*mean(projects) + ...
    0.15*exam1 + 0.15*exam2;
  gradeneeded = (target-current)/0.2;
  if gradeneeded > 100
    printf("I need a %.3f on the final.\n", gradeneeded)
    printf("RIP\n\n")
  elseif gradeneeded < 0
    printf("I need a %.3f on the final.\n", gradeneeded)
    printf("that means you I am above my target...\n\n")
  else
    printf("I need a %.3f on the final.\n\n", gradeneeded)
  end
end
% test using your data
quizzes = [60,80,70,90,100,100,40,90];
projects = [96,85,90,100,48,100,80,89,100];
exam1 = 89;
exam2 = 68;
% I want an 80
target = 80;
what_do_i_need_on_my_final(target, quizzes, projects, exam1, exam2)
final = 79;
[percent, letter] = wHaTisMygRAdE(quizzes, projects, exam1, exam2, final);
printf("My grade is %.4f (%s)\n", percent, letter)
printf("\n\nVectorization\n")
printf("=============\n\n")
% say we want to apply sin() function to all these x data points below
x = 0 : pi/2 : 2*pi;
% we can allocate our memory for the result vector and use a for loop
%  length() function get how many elements are in the array
y = zeros(1,length(x));
for i = 1:length(x)
  y(i) = sin(x(i));
end
y  % showing result to command window
% or we can simply do
z = sin(x);
z
% Octave recognizes x is an array and does this for loop in the background
% but actually in a much more efficient (parallel) way that is out of the scope
% of this class to explain.
% other vectorized examples, simple if doing operations with scalars
y = z + 100
y = z * 100
% operations with arrays and arrays have caveats
% back to our dot product example
x = [1; 2; 3];  % row vector
y = [4, 5, 6];  % column vector
w = y*x  % dot product gives us a scalar
% what if we simply wanted to multiply each corresponding element in the arrays?
% such as w = [y(1)*x(1), y(2)*x(2), ... ] to get back another array
% both arrays must be of the same dimension (either row vector or column vector)
% so we can transpose (switch the row and columns) of one of them
x
x = transpose(x);
x
% now do the element-wise computation
w = y .* x
% this gives us a matrix because 
% put a . in front of the operator, this is called an element-wise operation
printf("\n\nPlotting\n")
printf("========\n\n")
printf("DON'T LOOK HERE, THE PLOTS OPEN IN A SEPARATE WINDOW.\n\n")
% create data
t = 0:0.1:4;
hope = cos(t*pi/8);
% simply call the plot(x,y) function
plot(t, hope)
title("Your spirit over the course of an engineering degree.")
xlabel("years")
ylabel("spirit")
diary off  % stop saving output to file
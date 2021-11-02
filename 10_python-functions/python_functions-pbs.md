# Python Functions

## Define a function

A function is defined with *def*, followed by the name of the function and parenthesis.

As usual, we then indent anything we want to code as part of the function we have just defined.

The command *return* specifies the value we want the function to produce and pass outside the function.

```
def agree():
	return True
```

We can call a function, by typing its name followed by the parenthesis:

```
agree()
```

## Arguments


### Basics 

Values that are given *to* the function are called *arguments*. When you call a function with arguments, the values of those arguments are copied to corresponding *parameters* inside the function.

```
def echo(text):
	print("I have received: " + text)
```

Let's call it:

```
echo('random')
```

### Positional or Keyword arguments

Arguments are handled in a very flexible way in Python: valued passed to a function are copied to the corresponding parameters *in order*.

However, one can also pass the arguments by name: this way one avoids confusion with positions.

Let's see an example:

```
def write_menu(entree, main, dessert):
	print("\nToday's menu is going to be:")
	print("------------------------------")
	print("For starters, we're gonna have: " + entree + "\n---")
	print("As main, the Chef has chosen: " + main + "\n---")
	print("And finally for pudding, you'll enjoy: " + dessert + "\n======")
```
We can call the function by:

```
write_menu('crostini', 'bistecca', 'crostata')
```

But clearly if you write:

```
write_menu('crostata', 'bistecca', 'crostini')
```
You are not going to get the menu you'd ideally expect.


So one can solve this by using the keyword arguments:

```
write_menu(dessert='crostata', main='bistecca', entree='crostini')
```

We don't need to write the arguments in order anymore, as they appear in the function definition, because we have passed their value according to the arguments keywords.

However, if we do not pass all the arguments a function expects, we will receive an error:

```
write_menu(dessert='crostata', main='bistecca')
```

Using keyword arguments also allows us to define *default values*, which could be important.

```
def write_menu(entree='crostini', main='insalata', dessert='tiramisu'):
	print("\nToday's menu is going to be:")
	print("------------------------------")
	print("For starters, we're gonna have: " + entree + "\n---")
	print("As main, the Chef has chosen: " + main + "\n---")
	print("And finally for pudding, you'll enjoy: " + dessert + "\n======")
```

Now we don't need to pass all the arguments anymore.

```
write_menu(dessert='crostata', main='bistecca')
```


# Writing your code

## Script structure

Python is a very flexible language, and you don't need to stick to a specific structure to write functional code.

However, following a suggested structure helps preventing a few issues, especially when you'd like to re-use your code, and import the functions you have already written in order scripts.

### Begin with a shebang

A shebang is a commented line of code, which tells any interpreter what language is going to be present in the script you are writing.
This is the very first line of your script, in all languages you might be writing (of course, change the language if not Python):

```
#!/usr/bin/env python3
```

### Import any package you need

Following the shebang, we might need to import some packages useful for our scripts.
In this specific case, if we want to print the menu, we need to capture use inputs. 

Scripts have arguments, like functions do.
Capturing keyword arguments in a script is requires a package called *argparse* and a bit of learning, so to begin with we will stick with the *sys* package, which allows us to capture just arguments by position.

```
import sys
```


### Define your functions

If we wanted to create a script to write the menu, we would continue the script with the function we have defined:


```
def write_menu(entree='crostini', main='insalata', dessert='tiramisu'):
	print("\nToday's menu is going to be:")
	print("------------------------------")
	print("For starters, we're gonna have: " + entree + "\n---")
	print("As main, the Chef has chosen: " + main + "\n---")
	print("And finally for pudding, you'll enjoy: " + dessert + "\n======")
```

We can also define a function to capture the arguments we pass when we execute the script, assuming that by position they are going to be the entree, main and dessert.

```
def get_arguments():
	entree = sys.argv[1]
	main = sys.argv[2]
	dessert = sys.argv[3]
	return (entree, main, dessert)
```


### Define the main code

This is the important bit: once you have defined all the functions, there's a special *function* you should call *main()* (NB: no arguments), which defines what your script is going to do.

In our case, just call the menu function:

```
def main():
	entree, main, dessert = get_arguments()
	write_menu(entree, main, dessert)
```

### Setup the code execution

This is a special way to execute the function we have defined as *main()*, which most beginners do not use.

It is important however to acquire this practice early on, because it ensures you can re-use all your scripts, and that when you use the *import* statement on your code, it does **not** get executed:

```
if __name__ == '__main__':
	main()
```

Variables called with double underscores are special variables, and the code above is simply saying: if the name of the script I'm writing on the terminal corresponds to the name of the script file, then execute the main function. It basically detects when the code is executed, as opposed to when the code is imported.

Let's put all together in a file, called *menu.py* and let's execute it. An example below:

```
python menu.py zuppa stinco budino
```


# Problem solving & exercises

Now, in order to put all of this in practice, let's solve the first Rosalind Challenge.

## Count base frequency

The problem is the following:

- your script receives a DNA sequence
- your script results in the count of each base, i.e. how many A, T or C or G the sequenced you passed has.



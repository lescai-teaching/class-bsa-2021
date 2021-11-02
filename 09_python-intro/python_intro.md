# Introduzione Python


## Variables

Like in any other language, variables are *names* for values stored in your computer memory.

Variables are case-sensitive, i.e. *thing* is different from *Thing* and they're both different from *THING*.
Variable names cannot begin with a digit.


## Assignment

In Python you use *=* to assign a value to a variable

```
num = 2
name = "Joe"
```

## Data Types

There are several different data types, let's see a few of them

### Boolean

This is the typical *true* or *false*

```
test1 = True
test2 = 'true'
type(test1)
type(test2)
```

The above shows that the correct boolean value is with a capital letter. When I write with lowercase, it is not a boolean and it is going to be considered a string: for this reason I have to use the quotes.

If you type

```
test1 = true
```
it is going to throw an error, because it is looking for a variable name.

### Numbers

There is **no** *number* data type: the title of this paragraph serves to indicate several data types containing digits, used for arithmetics.

- integer (*int*): refers to numbers without decimal digits
- floating point (*float*): refers to numbers with decimal digits or scientific exponentials
- complex: refers to complex numbers

```
a = 2
b = 5
```


### String

Strings are sequences of characters, and they are created by using quotes.

```
data = 'this is a sentence'
data2 = "this is a sentence"
print(data)
print(data2)
```

There are single quotes ' and double quotes ". As you can see they are roughly equivalente, but they can be used together in order to create strings containing quotes.

```
sentence = "'awesome work' said mr. Ripley"
print(sentence)
```

Quotes enclosed in other quotes are not considered special characters anymore.

A string class has several properties, which is best to learn when needed. Just a few examples:

```
sentence = "this is just an example"
sentence.capitalize()
```

This has capitalised only the first word. You can use the following to capitalise every word:

```
sentence.title()
```

or every letter:

```
sentence.upper()
```


### List

The concept of a list should be relatively simple. It is a sequence of elements, like in many other languages.
They are assigned by using square brackets [ ... ], and the elements are separated by a coma.

They can be lists of any types.

```
songs = [ 'songA', 'songB', 'songC' ]
numbers = [ 1, 2, 3, 4, 5]
print(songs)
print(numbers)
```

### Tuple

While *lists* can be modified, in Python there's a special type of collection of elements which *cannot* be modified: this is very useful when you need to create data safe to be used, with a guarantee their structure, properties and content will not change.

Tuples are assigned using round brackets ( ... ), and like lists the elements are separated by comas.

```
collection = ( 'cA', 'cB', 'cC' )
num_collect = ( 1, 2, 3, 4 )
print(collection)
print(num_collect)
```


### Dictionaries

Dictionaries are a special kind of list, where an element is identified by a key. We call this **key-value** pair.
In this way, I can call a key of a dictionary and get its value. 
Another important difference with a list is that dictionaries are not ordered.

They are built using curly brackets as below, and the key-value indicated by colon, while elements separated by coma:

```
dictionary = {
	'day': 'a period of light',
	'night': 'a period in absence of light'
}
```

You can then call a value of a dictionary by its key, with the format below:

```
print(dictionary['day'])
```

or using the method *get*

```
dictionary.get('day')
```

You can also list all the keys with:

```
dictionary.keys()
```

or all the values:

```
dictionary.values()
```


## Adding & concatenating

Objects of the same type can be combined just using a *+* sign.

```
sentenceA = "this is only "
sentenceB = "an example"
sentence_complete = sentenceA + sentenceB
print(sentence_complete)
```

Same goes with numbers

```
a = 5
b = 2
c = a + b
print(c)
```

or lists:

```
list_1 = ['A', 'B']
list_2 = ['c', 'd', 'e']
list_all = list_1 + list_2
print(list_all)

```

But you cannot combine different types:

```
test = sentenceA + a
```

Unless you convert the number to a string:

```
test = sentenceA + str(a)
print(test)
```


## Slicing

Slicing means to get only part of a group of values.
It mainly applies to lists, but strings are also considered lists of characters.

Let's see a few examples:

```
list_a = ['dog', 'cat', 'parrot', 'elephant', 'eagle', 'mouse', 'rat' ]
new_sentence = "this is a series of random words"
```

Slicing is done using square brackets, and with the following format:

[ *start* : *end* : *step* ]

any of the three can be omitted.

try the following:

```
list_a[1:3]
new_sentence[6:9]
```

From the example above, you will see a *very* important characteristic of Python: **the first element of a list has index 0**.
It is important you remember this, because in R the first element of a list has index 1.

```
list_a[0:2]
new_sentence[5:8]
```
The next thing you will notice, is that the *end* means an *end offset minus 1*, or in other words that the end is **not** included.


## Conditions

Conditions are elements of code where we check if something is *true* or *false*.

The code uses the statement *if* / *else* and performs something in between.

### Syntax and Indenting

By indenting code, we mean the number of spaces from the beginning of the row: python uses indenting to separate blocks of code that need to be run together (loops, conditions) or that are dependent on the lines above (methods).

It is important to highlight this, because *not* indenting will prevent your conditionals and loops to run as they should.

Let's start with an if/else example:

```
animal = "dog"

if (animal == "cat"):
	print("this is a cat")
else:
	print("this is NOT a cat")

```

## Loops

### For and in

Loops are elements of code executed at every element of a collection (a list, a string, a dictionary):

```
word = "animal"
animals = ['cat', 'dog', 'parrot']
```

Now we can loop either on a string:

```
for letter in word:
	print(letter)
```

or through the list:

```
for element in animals:
	print(element)
```

When you have a loop it is also important to skip one round or even stop the loop, depending on some conditions. 

- *break* stops the loop
- *continue* skips the round

let's see an example:

```
for letter in word:
	if letter == "n":
		continue
	print(letter)
```

or if I want to exit the loop:

```
for letter in word:
	if letter == "m":
		break
	print(letter)
```


# Esercizi

Cerchiamo di riepilogare quanto visto fino ad ora.

Scrivi un codice che:
- crei una lista con le seguenti parole: uno, due, tre
- esegua un loop di questa lista
- all'interno del loop controlli se la parola corrisponde a *due* e nel caso scriva *si*, altrimenti scriva *no*



# Additional

If there is still time at the end of the class:

- slicing: omitting start/end
- slicing: step
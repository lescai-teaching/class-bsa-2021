# Biopython: working with NGS reads

From the entire biopython library we start by importing the methods to handle sequences, called **SeqIO**


```
from Bio import SeqIO
```

## Importing the reads

We have an example FASTQ in our repository: we use the method *parse* to read it:

```
data = SeqIO.parse( 'python_exercises/sample_seq.fastq', 'fastq')
```

Like for GenBank, we now have an iterable object which is a little difficult to inspect. We can look at the methods to access it by typing:

```
var(data)
```

Since it is not a list, but a complex iterable object, we can just pop out one element with *next* and have a look at it:

```
record = next(data)
```
now we can print it:

```
print(str(record))
```

Or print more details:

```
print(record.id, "\n", record.description, "\n", record.seq)
print(str(len(record.seq)))
```

The code below is ok because we don't have many sequences, but it will fill in your page and it is not recommended for real data:

```
for record in data:
    print(str(len(record.seq)))
```



## Counting the reads

We can use the same loop in order to count the reads with a counter.

```
count = 0
for rec in SeqIO.parse('python_exercises/sample_seq.fastq', "fastq"):
    count += 1
print("%i reads" % count)
```

Obviously, there are smarter ways to do that.


## Filtering the reads

We can access the Phred quality of the reads and use it to filter our data: this is an important real case.

```
good_reads = (rec for rec in \
              SeqIO.parse('python_exercises/sample_seq.fastq', "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 10)
count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
print("Saved %i reads" % count)
```

These are simulated reads, and their quality is clearly not great.
Let's get some more details.

## Inspecting Quality

We are doing a bit of math here, so we need to import the appropriate libraries.
```
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pyplot
import numpy as np
```

In order to save these data for future use, we instantiate a few empty lists which we can fill with loops:
```
pList = []
phredList =[]
lengths = []
```


It is important to notice that we need to parse the fastq file again because the list is consumed every time we loop through it.

```
data = SeqIO.parse( 'python_exercises/sample_seq.fastq', 'fastq')
```

Now we can loop through it and extract the data.
```
for record in data:
    pArray = []
    phredArray = []
    qualities = record.letter_annotations["phred_quality"]
    for Q in qualities:
        # convert the PHRED score to a probability
        p = 10**(-float(Q)/10)
        # append this specific probabity to array for this sequence
        pArray = pArray + [p] 
        phredArray += [Q]
    # now append the sequence's probablity to a list of all of them    
    pList = pList + [pArray]
    # also append the phred qualities as they are
    phredList += [phredArray]
    # store the length of the probability array (same as length of sequence)
    lengths += [len(pArray)]
```

we can see lenghts are not all the same

```
min(lengths)
max(lengths)
```

now for plotting. Initialize x and y arrays to plot:

```
x = []
y = []
```

add the average probabilties to the y values:

```
for i in range(min(lengths)):
    p_i = []
    for p in pList:
        p_i = p_i + [p[i]]
    pAvg = sum(p_i)/len(p_i)
    x.append(i)
    y.append(pAvg)
```

plot the x and y values:

```
pyplot.figure(0)
pyplot.plot(x,y)
pyplot.xlabel('position (nt)')
pyplot.ylabel('Averge error probability')
pyplot.savefig('quality.png')
```


# Pandas

Pandas is an important tool for data science in python, because it allows to store data in a tabular format.

This is much easier to read and understand than the array structure or the list of lists structure we have used above.


```
import pandas as pd
```

## Convert array into data frame

Now we can use Pandas to convert our list of lists to a data frame.

```
df = pd.DataFrame(pList)
```
If we print it, we have a more familiar data structure in front of us.

```
print(df)
```

We can use the method *shape* to see the number of rows and columns of our table.

```
df.shape
```

Pandas Data frames are like dictionaries, where columns could have names.

In our case we did not assign names, so each column can be accessed by number like below:


```
df[0]
```

We will do the same for plotting average probabilities, but this time with much less code:

```
base = []
val = []
for i in range(min(lengths)):
    base.append(i)
    val.append(sum(df[i])/len(df[i]))
```

Now we can plot the base position in x and average p in y values:

```
pyplot.figure(1)
pyplot.plot(base,val)
pyplot.xlabel('position (nt)')
pyplot.ylabel('Averge error probability')
pyplot.savefig('quality_2.png')
```


We could have done the same for the phred qualities

```
phredData = pd.DataFrame(phredList)

print(phredData)
```

Like above we instantiate our list of values to be plotted and we fill it in, with a loop

```
base = []
val = []

for i in range(min(lengths)):
    base.append(i)
    val.append(sum(phredData[i])/len(phredData[i]))
```

Now we are ready to plot the base position in x and average p in y values:

```
pyplot.figure(2)
pyplot.plot(base,val)
pyplot.xlabel('position (nt)')
pyplot.ylabel('Averge Phred Quality')
pyplot.savefig('quality_phred.png')
```


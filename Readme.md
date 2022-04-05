# Lineage
*Ancestral recombination graph in Asymptote made easy*

Usable, but still a work in progress.

## Installation
Simply add `Lineage.asy` to your favorite [Asymptote search
path](https://asymptote.sourceforge.io/doc/Search-paths.html). The library can
be imported with a simple `import Lineage` statement.

## Example
This is an example of a simple ancestral recombination graph outputed by Lineage:

![Lineage example](/.readme_assets/example.svg)

It has been generated by the following code:
```asymptote
import Lineage;

settings.outformat = "svg";

unitsize(1cm);

bool[] type1 = {true, false, false, false, false};
bool[] type2 = {true, false, true, false, true};
bool[] type3 = {true, false, true, false, false};
bool[] type4 = {false, false, false, false, false};

Arg arg = Arg(type1);
arg.newleaf(type2, 2);
arg.newleaf(type3);
arg.newleaf(type4, 2);

arg.coalesce(1, 2);
arg.coalesce(3, 6);
arg.coalesce(4, 5);
arg.recombine(7, 1);
arg.coalesce(8, 10);
arg.coalesce(0, 9);
arg.coalesce(11, 12);

arg.draw();
```

## TODO
- [ ] Write documentation
- [ ] Write tutorial
- [ ] Make vertices latitude customizable
- [ ] Allow multiple mutations per coalescence event
- [ ] ...

## Acknowledgements
- [Asymptote](https://asymptote.sourceforge.io/): The Vector Graphics Language.

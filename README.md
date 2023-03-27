# EVOLVE

A genetic algorithm tool to mutate protein structures based on Richardson's Rotamer library and to
fitness evaluations

## Citations

Please cite the following work if you use and modify EVOLVE to add your own evaluators. 

`EVOLVE citation`

## Getting Started

These instructions will get you a copy of the project up and running on your local machine. 

### Prerequisites

EVOLVE is python3 based and depends on the following packages

```
pygmo
amber18
openmm
openbabel
```

### Installing

A step by step series of examples that tell you how to get a development env running

tbc

## Documentation 

Please find detailed instructions on the different functions of the code, how to implement a new
evaluator and a detailed instructions to run the tests in the [Documentation](http://lcbc.epfl.ch)

## Built With

* [openbabel](http://openbabel.org/) - Openbabel
* [Amber](https://ambermd.org/) - to build toplogies
* [OpenMM](http://openmm.org) - to run molecular dynamics

## Contributing

Please read `CONTRIBUTING.md` for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Nicholas J. Browning** - *code development, mutation modules, stability indicator*
* **Justin Villard** - *amber gpgg* 
* **Marta Gomez** - *Initial thermostability code* -
* **Simon L. Duerr** - *Multi-objective optimization routines, tools for enyzme design, code development*
* **Ursula Roethlisberger** - *Prinicipal investigator*


See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details

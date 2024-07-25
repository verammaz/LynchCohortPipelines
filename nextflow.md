# nf-core/sarek 

Visit the official site [here](https://nf-co.re/sarek/3.4.2/) for full details about the nextflow nf-core/sarek pipeline

## Usage on Minerva

1. Navigate to the directory where you want your nextflow executable to be located.
2. Run the following:
```bash
wget -qO- https://get.nextflow.io | bash
```
3. Check the version:
```bash
/path/to/nextflow -v
```
4. Edit the file ~/.nextflow/assets/nf-core/sarek/nextflow.config by adding the following chunk:
```bash
executor {
    name   = 'local'
    cpus   = 48
    memory = '64GB'
}
```


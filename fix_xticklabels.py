path = '/scratch/easmit31/mitotyping/mitotyping_plots.py'
txt  = open(path).read()

old = '''    mat.columns = [f"{str(a)}\\n({float(age):.0f}y)"
                   for a, age in zip(animal_means['animal_id'], animal_means['age'])]'''

new = '''    mat.columns = [f"{str(a)} ({float(age):.0f}y)"
                   for a, age in zip(animal_means['animal_id'], animal_means['age'])]'''

if old in txt:
    txt = txt.replace(old, new)
    open(path, 'w').write(txt)
    print("Fixed x tick labels — now on one line")
else:
    print("Pattern not found — searching for similar lines:")
    for i, line in enumerate(txt.split('\n')):
        if 'animal_id' in line and 'age' in line and 'columns' in line:
            print(f"  {i}: {line}")

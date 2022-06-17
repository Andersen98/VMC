# Example Scheme

Remark: We will just link to the pyscf site so we can reference their sphinx docs. Following the guide [from the sphinx docs.](https://docs.readthedocs.io/en/stable/guides/intersphinx.html)

Put this in the `conf.py` file
``` python
# conf.py file

extensions = [
    'sphinx.ext.intersphinx',
]


intersphinx_mapping = {
    'pyscf': ('https://pyscf.org/', None),
}
```

Now we can use the `sphinx` name with a cross-reference rule 
```rst
- :ref:`sphinx:ref-role`
- :ref:`:ref: role <sphinx:ref-role>`
- :doc:`sphinx:usage/extensions/intersphinx`
- :doc:`Intersphinx <sphinx:usage/extensions/intersphinx>`
```
```markdown
- {ref}`sphinx:ref-role`
- {ref}`:ref: role <sphinx:ref-role>`
- {doc}`sphinx:usage/extensions/intersphinx`
- {doc}`Intersphinx <sphinx:usage/extensions/intersphinx>`

```
- `${method_name}`
  - `README.md`
  - `prepVMC.py`
  - `${example_system_1}`
    - `${m}z`
      - `${example_system_1}.py`
      - `dice.dat`
      - `dice.out`
      - `dqmc.json`
      - `dqmc.out`
      - `pyscf.out`
    -`${m+1}z`
        `...`
  - `${example_system_2}`
    - `...`
  - `...`
  - `${example_system_N}`
    - `...`


We describe the elements in the scheme in turn below.

- `${method_name}`: name of method used for the examples. 
- `${README.md}`: lists the software used and their respective sources.
- `${example_system_x}`: a folder that holds the example for a system. Could be benzene for example. 
- `${prepVMC.py}`: A utility file for prepping the example. You import it as a module using `import prepVMC`. The module contains the following functions:
    - `prepVMC.doRHF(**kwargs)`, `prepVMC.doUHF(**kwargs)`, `prepVMC.doGHF(**kwards)`,`prepVMC.localizeAllElectron(**kwargs)`, `prepVMC.localizeValence(**kwargs)`, `bestDetValence(**kwargs)`, `prepVMC.writeFCIDUMP(**kwargs)`, `basisChange(**kwargs)`, `writeMat(**kwargs)`, `readMat(**kwargs)`,`makeAGPFromRHF(**kwargs)`, `makePfaffFromGHF(**kwargs)`,`addNoise()` , `prepAllElectron()`,`prepValence`, `write_dqmc`, `write_ccsd`, `write_uccsd`, `findSiteInUnitCell`,`findSiteAtRowNCol`, `findNeighbors`,  - Where the argument is a {ref}`pyscf:module-pyscf.gto.molek`.
    
```{eval-rst}
.. py:function:: prepVMC.doRHF(mol)

    Given a molecule, does RHF and returns the result.
    
    :param mol: Molecule to do RHF on.
    :tpye mol: :ref: `pyscf:module-pyscf.molek`
    :rtype: module-pyscf.scf.hf
```

## To do
- [ ] Re-write constructors to handle data structures for ArupCompute
- [ ] Add error handling for XTRACT
        - Clean reporting and termination if XTRACT is not found on system
        - Write temporary data to appropriate Windows directory
- [ ] Fix working directory limitations in hysBeams
- [ ] Add option to start with .xpj files


## Validations
- [x] Check XTRACT outputs implementing Mander versus not 
- [x] Double check class(info) values are not erroneously instantiated in the middle of analyses
- [x] Fix EC2 axial calculations to remove potential for unreasonable shear strains


## Wish list
- [ ] Remove XTRACT dependency
- [ ] Add option for L-beams
- [ ] Simplify inputs


## Done
- [x] Setup iterative solution for envelope of P-M interaction based on a broad range of yield strains

        - Apply to 'Unconfined' and 'Confined' compression limiting strains
        - Iterate from 1.2E-3 to xxE-3 
            + based on EC2 
            + user inputs nominal f'c from drawings
            + user defined upper bound to be replaced in later stages by automatic calculation

- [x] Replace unnecessary classes with methods and restructure pipeline accordingly
- [x] Add method for computing moment-plastic rotation according to ASCE 41 (e.g. Tables 10-7, 10-8, 10-9; @Yuli.Huang)

        - Parse M-Phi bilinearization from XTRACT
        - Determine failure strain according to code
        - Update card writing module to write the appropriate values to PR1, PR2, PR3, or PR4 (choice based on consequence class)

- [x] Introduce git/config settings to handle user-specific settings
- [x] Build an auto-update module
- [x] Add method to write FEMA rotation thresholds based on ASCE 41-17 to material cards
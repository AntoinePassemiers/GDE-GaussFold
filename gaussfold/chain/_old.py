    def create_all_atom_model(self, chain, gds, ssp, contact_threshold):
        L = len(ssp)
        segment_ids = [0]
        for i in range(1, L):
            if ssp[i] != ssp[i-1]:
                segment_ids.append(segment_ids[i-1] + 1)
            else:
                segment_ids.append(segment_ids[i-1])

        atoms = chain.atoms()
        model = AllAtomModel(chain)

        for i in range(len(chain)):
            for j in range(max(0, i-self.sep)):
                if gds[i, j] == 1: # Contact
                    #a = (chain[i].CB if not isinstance(chain[i], Glycine) else chain[i].CA).identifier
                    #b = (chain[j].CB if not isinstance(chain[j], Glycine) else chain[j].CA).identifier
                    a = chain[i].CA.identifier
                    b = chain[j].CA.identifier
                    model.add_distance_restraint(a, b, lb=3.5, ub=8.) # TODO

        for bond in chain.bonds():
            i = bond.atom1.identifier
            j = bond.atom2.identifier
            mu = bond.distance
            print(mu)
            model.add_distance_restraint(i, j, lb=mu-0.1, ub=mu+0.1)

        for i in range(L):
            if ssp[i] == 0: # Helix
                model.add_angle_restraint(i, psi_lb=-46.4, psi_ub=-36.6, phi_lb=-68.1, phi_ub=-58.9)
            elif ssp[i] == 1: # Strand
                # TODO: make distinction between parallel and anti-parallel beta-sheets
                model.add_angle_restraint(i, psi_lb=122.6, psi_ub=145.6, phi_lb=-131.9, phi_ub=109.9)

        return model
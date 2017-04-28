# Index parameters
indparams = {'NaID_classic':
                 {'w1': {'ww0': 5895.92, 'dww': 0.5},  # Core D1
                  'w2': {'ww0': 5889.95, 'dww': 0.5},  # Core D2
                  'wr1': {'ww0': 5805.0, 'dww': 10.0},
                  'wr2': {'ww0': 6090.0, 'dww': 20.0}
                  },

             'NaID_telluric':
                 {'w1': {'ww0': 5895.92, 'dww': 0.5},  # Core D1
                  'w2': {'ww0': 5889.95, 'dww': 0.5},  # Core D2
                  'wr1': {'ww0': 5805.0, 'dww': 0.0},
                  'wr2': {'ww0': 6090.0, 'dww': 20.0}
                  },

             'Halpha_narrow':
                 {'w1': {'ww0': 6562.808, 'dww': 0.6},  # Boisse+2009
                  'wr1': {'ww0': 6605.0, 'dww': 20.0},  # Cincunegui+2007
                  'wr2': {'ww0': 6605.0, 'dww': 0.0}
                  },

             'Halpha_broad':
                 {'w1': {'ww0': 6562.808, 'dww': 1.6},  # DaSilva+2011
                  'wr1': {'ww0': 6605.0, 'dww': 20.0},  # Cincunegui+2007
                  'wr2': {'ww0': 6605.0, 'dww': 0.0}
                  },

             'Halpha_boisse':
                 {'w1': {'ww0': 6562.808, 'dww': 0.68},  # Boisse+2009
                  'wr1': {'ww0': 6550.0, 'dww': 10.76},  # Boisse+Girault+Rey
                  'wr2': {'ww0': 6580.0, 'dww': 8.75}
                  },

             'Halpha_dasilva':
                 {'w1': {'ww0': 6562.808, 'dww': 1.6},
                  'wr1': {'ww0': 6550.87, 'dww': 10.75},
                  'wr2': {'ww0': 6580.31, 'dww': 8.75}
                  },

             'CaII_wilson':
                 {'w1': {'ww0': 3933.664, 'dww': 2.0},
                  'w2': {'ww0': 3968.470, 'dww': 2.0},
                  'wr1': {'ww0': 4001.07, 'dww': 20.0},
                  'wr2': {'ww0': 3901.07, 'dww': 20.0}
                  },

             'CaII_lovis':
                 {'w1': {'ww0': 3933.664, 'dww': 2.0},
                  'w2': {'ww0': 3968.470, 'dww': 2.0},
                  'wr1a': {'ww0': 3918.664, 'dww': 20},
                  'wr1b': {'ww0': 3948.664, 'dww': 20},
                  'wr2a': {'ww0': 3953.470, 'dww': 20},
                  'wr2b': {'ww0': 3983.470, 'dww': 20},
                  }
             }
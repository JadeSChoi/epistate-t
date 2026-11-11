"""Minimal smoke check of scope, missing dates, and the effect of the example's reader branch."""
from demo import answer,records,Config,simulate
def test():
    rows=records([.2,.55,.21],'constructed-culture','made-up-site','treated')
    assert answer(rows,'series_return',origin='simulation')['value'] is True
    assert answer(rows,'cell_return',origin='simulation')['reason']=='not_tracked'
    assert answer(rows,'series_rate',origin='simulation')['reason']=='no_duration'
    assert answer(rows,'series_two_visit_return',origin='simulation')['reason']=='too_few_visits'
    off=simulate(Config(seed=20261106,retention=1,read_strength=0),detail=False)['rmse']
    on=simulate(Config(seed=20261106,retention=1,read_strength=4),detail=False)['rmse']
    assert (off[:3]==on[:3]).all() and on[-1]<off[-1]
    print('Scope, missing time, visit count and illustrative simulation checks passed.')
if __name__=='__main__': test()

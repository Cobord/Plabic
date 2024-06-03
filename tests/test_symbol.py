"""
test for symbol alphabet
"""

from plabic.symbol_alphabet import Symbol,SymbolWord

def test_0log():
    """
    zero tensorands
    """
    errored = False
    try:
        _nonexistent = SymbolWord([],lambda z:[z],lambda z: {z:f"1/{z}"})
    except ValueError:
        errored = True
    assert errored

def test_log():
    """
    symbol is just a single letter for a single log
    """
    loga = SymbolWord(["a"],lambda z:[z],lambda z: {z:f"1/{z}"})
    assert loga.my_variables() == frozenset(["a"])
    da_over_a = loga.to_integrands()
    assert da_over_a == [{"a":"1/a"}]
    _five_loga = Symbol([(5,loga)])

    logalpha = SymbolWord(["alpha"],lambda z:[z],lambda z: {z:f"1/{z}"})
    assert logalpha.my_variables() == frozenset(["alpha"])
    dalpha_over_alpha = logalpha.to_integrands()
    assert dalpha_over_alpha == [{"alpha":"1/alpha"}]
    _eight_logalpha = Symbol([(8,logalpha)])


def test_2log():
    """
    two tensorands
    int dalpha/alpha int db/b
    log alpha log b + C log alpha + D log b + CD
    """
    logalphab = SymbolWord(["alpha","b"],lambda z:[z],lambda z: {z:f"1/{z}"})
    assert logalphab.my_variables() == frozenset(["alpha","b"])
    diff_forms = logalphab.to_integrands()
    assert diff_forms == [{"alpha":"1/alpha"},{"b":"1/b"}]

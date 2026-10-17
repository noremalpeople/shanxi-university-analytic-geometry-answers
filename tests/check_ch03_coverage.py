"""Prevent accidental omissions of chapter-three numbered exercises."""
import re
from pathlib import Path

root = Path(__file__).resolve().parents[1] / 'chapters' / '03' / 'chapters'
expected = {'03-01.tex': 16, '03-02.tex': 11, '03-03.tex': 7,
            '03-04.tex': 24, '03-05.tex': 5}
total = 0
for filename, maximum in expected.items():
    tex = (root / filename).read_text(encoding='utf-8')
    actual = list(map(int, re.findall(r'\\Q\{(\d+)\}', tex)))
    required = list(range(1, maximum + 1))
    assert actual == required, f'{filename}: expected {required}, got {actual}'
    assert tex.count(r'\ans{') >= maximum - 1, f'{filename}: missing conclusion markers'
    total += len(actual)
assert total == 63
print('PASS: all 63 Chapter 3 numbered exercises are covered in source order')

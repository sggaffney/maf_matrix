__author__ = 'sgg'

from unittest import TestCase

# import maf_matrix


class Test(TestCase):
    def test_is_string(self):
        s = 'Remember to write test'
        self.assertTrue(isinstance(s, basestring))

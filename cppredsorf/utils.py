import pkg_resources


def get_model_range_hexamer():
    '''
    '''

    m_f = pkg_resources.resource_filename(__name__, 'Model/Integrated.model')
    r_f = pkg_resources.resource_filename(__name__, 'Model/Integrated.range')
    h_f = pkg_resources.resource_filename(__name__, 'Model/Intwgrated_Hexamer.tsv')
    return m_f, r_f, h_f

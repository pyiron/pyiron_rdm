def get_dataset_assignments(o):
    return {i.code: i.get_property_assignments().df for i in o.get_dataset_types()}


def get_object_assignments(o):
    return {i.code: i.get_property_assignments().df for i in o.get_object_types()}


def get_controlled_vocabulary(o):
    return {i.code: i.get_terms().df for i in o.get_vocabularies()}

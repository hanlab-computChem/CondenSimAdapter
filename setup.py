from setuptools import setup, find_packages

setup(
    name='CondenSimAdapter',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'CondenSimAdapter': ['**/*'],
    },
    description='An adapter for CG to allatom protein condensate simulation.',
    author='Xiaojing Tian',
    author_email='tianxj15@tsinghua.org',
    entry_points={
        'console_scripts': [
            'adapter=CondenSimAdapter.cli:main',
        ],
    },
    install_requires=[
        'click>=8.0',
        'pyyaml>=6.0',
    ],
)


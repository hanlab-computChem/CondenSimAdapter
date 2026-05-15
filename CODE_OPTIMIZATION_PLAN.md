# CondenSimAdapter 代码工程质量优化计划

## Summary

目标是先把仓库从“能跑但维护成本偏高”优化到“静态检查可控、模块边界清晰、关键路径可安全修改”。本计划不优先做性能算法优化，除非重构过程中发现明显低风险问题。

当前发现：

- Python 包约 `22,407` 行，测试约 `5,051` 行。
- `ruff check CondenSimAdapter tests` 当前报 `1327` 个问题，其中 `997` 个可自动修复。
- 主要维护风险集中在 `minimize/softcore.py`、`minimize/topology_builder.py`、`forcefield/registry.py`、`core/entanglement.py`、`backmap/cg2all/*`。
- `pyproject.toml` 和 `pytest.ini` 都配置了 pytest，且 marker/addopts 不完全一致，需要合并来源。

## Public API / Interface Changes

- 保持现有公开 CLI 不变：`adapter cg/backmap/minimize/to_run/init/info/forcefield/models` 命令名、参数、输出目录结构不做破坏性修改。
- 保持现有 Python import 兼容：`CondenSimAdapter.cli.commands` 作为 backward compatibility 模块继续保留。
- 新增内部模块可以使用下划线命名或放在已有子包内，但不新增用户必须直接调用的新 API。
- 若拆分大文件，只通过原文件 re-export 或内部导入保持旧路径可用，特别是 `minimize.softcore` 和 `minimize.topology_builder`。

## Implementation Changes

### 1. 质量基线与配置统一

- 合并 pytest 配置：以 `pyproject.toml` 为唯一来源，迁移 `pytest.ini` 中的 `-W ignore` 和 marker，删除或停用重复配置。
- 明确 ruff 分层策略：
  - 项目自有代码逐步达到 `ruff check` 干净。
  - `CondenSimAdapter/backmap/cg2all/` 视为移植/第三方风格代码，先在 `pyproject.toml` 加 `per-file-ignores`，只修影响运行的 `F` 类错误，不做大规模风格重排。
  - `tests/` 自动清理 import、空白、未使用变量。
- 更新 `Makefile`：
  - `lint` 使用 `python -m ruff check` 和 `python -m ruff format --check`，避免环境 PATH 差异。
  - 新增 `test-fast`，运行不依赖 GPU/GROMACS/慢集成的单元测试子集。

### 2. 逐文件优化清单

- `CondenSimAdapter/minimize/softcore.py`
  - 保持 `GromacsTopFileWithSoftcore` 对外路径不变。
  - 拆出 3 个内部模块：GROMACS topology 解析、softcore force 构造、GBSA/CustomNonbonded 助手。
  - 给复杂表达式构造函数加最小单元测试，覆盖 force expression 字符串、lambda 参数、cutoff/switch 参数。
  - 不改变数值公式和默认参数。

- `CondenSimAdapter/minimize/topology_builder.py`
  - 拆分外部命令执行、序列/PDB/FASTA 解析、topology merge/verify 三类逻辑。
  - 将 `subprocess.run` 包装成单一 helper，统一 `capture_output=True`、错误信息和命令缺失提示。
  - 给 `_build_pdb2gmx_input`、HIS 计数、`verify_files_valid`、topology merge 增加或补齐 fixture 测试。
  - 保持现有 `run_pdb2gmx_for_topology`、`run_pdb2gmx_for_structure` 函数名可导入。

- `CondenSimAdapter/forcefield/registry.py`
  - 拆分 registry 数据模型、内置 force field 列表、自定义 force field 持久化。
  - 替换静默 `except Exception: pass` 为带上下文的 warning 或明确降级返回。
  - 对 JSON index 读写使用临时文件 + 原子替换，避免用户 forcefield registry 写坏。
  - 补齐测试：损坏 JSON、重复注册、删除不存在项、内置 force field 不可删除、自定义路径缺文件。

- `CondenSimAdapter/core/entanglement.py`
  - 保留现有几何算法结果，先做命名、类型、边界条件整理。
  - 将低层几何函数和报告聚合逻辑分开，避免 analyzer 直接承载所有细节。
  - 补齐边界测试：共面、端点相交、零长度 segment、周期边界 unwrap。

- `CondenSimAdapter/core/simulation.py`
  - 将 `_run_pipeline` 拆成 build、minimize、production、save、entanglement 五个私有步骤。
  - 将 progress/checkpoint 频率逻辑提成纯函数，测试 steps 小于、等于、大于 batch 的情况。
  - 保持 `CGSimulation.run()` 行为不变。

- `CondenSimAdapter/src/plumed_generator.py`
  - 统一 legacy/new component 字段解析，减少多处 `getattr` 和宽泛异常。
  - 将 fdomains 解析、global index map、contact map 生成分别测试。
  - 错误信息包含 component name、domain range、文件路径。

- `CondenSimAdapter/cli/__init__.py` 和 `CondenSimAdapter/cli/commands.py`
  - 保留 lazy command 机制。
  - 为 `_LazyCommand` / `_LazyGroup` 增加类型标注和轻量测试，确保 `adapter --help` 不导入 heavy deps。
  - `commands.py` 只作为兼容 re-export，不再承载新逻辑。

- `CondenSimAdapter/backmap/cg2all/*`
  - 不做大规模格式化，避免污染移植代码 diff。
  - 只修运行风险：未使用但会触发失败的变量、裸 `except` 中会吞掉关键文件损坏的路径、模型下载错误提示。
  - 对该目录设置 ruff 例外，后续若要深度维护再单独建专项。

- `tests/*`
  - 自动修复 import 顺序、空白行、未使用 import。
  - 将依赖 OpenMM/GROMACS/GPU 的测试 marker 明确化，避免普通 CI 误跑重依赖测试。
  - 把 smoke、unit、integration 命名边界整理清楚。

## Test Plan

- 基线记录：
  - `python -m ruff check CondenSimAdapter tests`
  - `python -m pytest tests/unit -ra`
- 每阶段验收：
  - hygiene 阶段：`python -m ruff check CondenSimAdapter tests` 问题数明显下降，且不引入行为变更。
  - refactor 阶段：每拆一个模块，运行对应单测，例如 `tests/unit/test_minimize_topology.py`、`tests/unit/test_core_entanglement.py`、`tests/test_forcefield_registry_cli.py`。
  - CLI 阶段：`python -m pytest tests/unit/test_cli_entrypoint.py tests/test_forcefield_registry_cli.py -ra`。
- 最终验收：
  - `python -m ruff check CondenSimAdapter tests`
  - `python -m ruff format --check CondenSimAdapter tests`
  - `python -m pytest tests/unit -ra`
  - 可选环境具备依赖时运行：`python -m pytest tests/integration -ra -m "not slow"`。

## Assumptions

- 优化重点采用“工程质量”，不是优先做运行性能 profiling。
- 粒度采用“逐文件清单”，但实施时按风险分阶段提交，避免一次性大 diff。
- `backmap/cg2all` 按移植代码处理，默认不做全量风格重写。
- 不破坏现有 CLI 和 Python import 路径；所有拆分都通过兼容导入保持旧调用可用。

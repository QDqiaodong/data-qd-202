# data-qd-202 · 晶体分数坐标与位点展开

验证显式有限仿射操作集合的闭包，展开规范分数坐标；同源特殊位置去重，不同源占位合并并核防超占位，输出晶胞组成和位点来源。

本仓主要交付为Python业务API与下游可消费的数据，不以通用CLI、报表或后台模板为产品。可观察字段、单位、身份、上限及失败边界见 [数据契约](docs/contract.md)。

## 独立安装与验证

Python >=3.11，实际验证CPython3.12.14、uv0.12.19。运行和测试仅标准库；构建固定setuptools80.9.0、wheel0.45.1。各仓独立 `.venv`，可以共享公共下载缓存，不共享site-packages。

从本仓根目录执行：

```sh
uv venv --offline --python 3.12 .venv
uv pip install --offline --python .venv/bin/python -r requirements-build.lock
uv pip install --offline --python .venv/bin/python --no-build-isolation --no-deps .
uv pip check --python .venv/bin/python
.venv/bin/python -I -m unittest discover -s tests -v
.venv/bin/python -I examples/run_sample.py
```

以上命令在本地现有公共缓存上真实执行通过，无需新增下载。新机器缓存为空时，先准备匹配Python与锁定公共构建包，或在可联网安装环境移除安装命令中的offline；API、测试与样例运行不需要网络。示例创建 `run/sample.json`，目标已存在会拒绝覆盖；重复运行请先自行保留产物并选择新的演示工作副本。

## 使用与实际结果

```python
from site_orbits import process
from site_orbits.io import load, publish

result = process(load("fixtures/input.json"))
publish("run/another-sample.json", result)
```

恒等和反演操作展开2个源位点得到3个位置，组成X=1、Y=1；原点特殊位点占位只计一次。

安装后的样例摘要：

```json
{"contributions": 3, "expanded_positions": 3, "operations": 2, "source_sites": 2}
```

已通过 6 项业务测试，另实际核对已有目标、符号链接目标、JSON重复键、非法业务输入和临时发布文件清理。数值和排序语义见契约，不以文件哈希代替业务结果断言。

## 当前边界与来源

不推断空间群，不读取真实材料库，不验证晶胞度量等距性或化学可行性，不做物理模拟。

全部输入样例原创合成，MIT许可，见 [来源](fixtures/SOURCES.md)。仅复用MIT通用构建、字段校验与原子新建IO脚手架；业务实现独立编写。当前精确公共依赖镜像准备查询结果为 `true`；这不影响已通过的宿主离线独立环境验证，也不构成容器验收记录。本次没有构建或启动容器。
